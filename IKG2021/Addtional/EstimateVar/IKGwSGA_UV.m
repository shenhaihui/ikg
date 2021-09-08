function regret_plot = IKGwSGA_UV(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2)

IKG_pra.SGA_m = SGA_m;       % batch size
IKG_pra.SGA_K = SGA_K;       % iteration number
IKG_pra.SGA_b1 = SGA_b1;     % b1 
IKG_pra.SGA_b2 = SGA_b2;     % b2; in SGA, b_k = b1/(k^(b2))

M = 5; % number of alternatives
% constraints of x
lb = zeros(d,1); % lower bound
ub = 10*ones(d,1); % upper bound

% theta, tau2 -  parameters for Gaussian Covariance (prior)
theta = 1 / d * ones(d,1,M);
tau2 = 1 * ones(1,1,M);

% pack fixed parameters for IKG
IKG_pra.theta = theta; IKG_pra.tau2 = tau2;
IKG_pra.lb = lb; IKG_pra.ub = ub;
IKG_pra.density_type = density_type;
IKG_pra.cost_type = cost_type;

L = 30; % replication of IKG
regret_cost_L = cell(1,L);
    
% Monte Carlo samples 
s_pseudo = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s_pseudo);

J2 = 1000*d^2; % for evaluating regret
X2 = rand_c(d,J2,density_type); 
J = 500*d^2; % for sample average approximation (SAA) of IKG

% random starts, not related to distribution of covariates
% generated out of the loop of sampling, for relatively fair comparision
v0_pool_L = rand(d,M,2*C,L)*10;

% ----- estimate sampling variance -----
var_tau2 = 1;  % parameters for Gaussian Covariance (prior) used in var est
var_theta = 1 / d * ones(d,1);
k = 6; % number of points to estimate sampling var
m = 10; % number of samples to estimate simulation var at each point
varPra.var_tau2 = var_tau2;
varPra.var_theta = var_theta;
varPra.k = k;
var_X_L = zeros(d,k,L);
var_Cov_Inv_y_M_L = zeros(k,M,L);
for ell = 1:L
    var_X = lhsdesign(k,d)*10; % LH design
    var_X = var_X';
    var_X_L(:,:,ell) = var_X;
    for i = 1:M
        noise = sqrt(noise_var(i,var_X));
        noise = repmat(noise,m,1);
        noise = noise .* randn(m,k);
        var_est = var(noise)'; % use sample variance as estimate
        CovM = zeros(k,k);
        for Xi = 1:k
            CovM(:,Xi) = cov_vector(var_X(:,Xi), var_X, var_theta, var_tau2);
        end
        var_Cov_Inv_y_M_L(:,i,ell) = CovM \ var_est;
    end
end
% var_estPra.var_X_L = var_X_L;
% var_estPra.var_Cov_Inv_y = var_Cov_Inv_y;
% --------------------------------------

% parellel computation
parpool(parpool_n);
parfor ell = 1:L 
% for ell = 1:L
    
var_X = var_X_L(:,:,ell);
var_Cov_Inv_y_M = var_Cov_Inv_y_M_L(:,:,ell);
    
s_pseudo = RandStream('mt19937ar','Seed',ell);
RandStream.setGlobalStream(s_pseudo);
    
regret_cost = zeros(2,2*C);
v0_pool = v0_pool_L(:,:,:,ell);
X = rand_c(d,J,density_type);

% to compute regret
theta_true = zeros(M,J2);
for i = 1:M
    theta_true(i,:) = griewank_rev(i,X2); 
end
mu_i_n_x0 = mu_prior(-M,X2)';
mu_i_n_x = mu_i_n_x0;
true_max = max(theta_true);
[~,chosen_index] = max(mu_i_n_x);
index = sub2ind(size(theta_true),chosen_index,1:J2);
true_chosen = theta_true(index);
regret_cost(1,1) = mean(true_max - true_chosen);

N0 = C*0.5;  % the actual used size for each alternative; to save memory.
V = zeros(d,N0,M); % visited points, for each slice (alternative), each colum is a point
yn = zeros(1,N0,M); % visited values, for each slice, each colum is an observation
mun = zeros(1,N0,M); % visited prior mean
Vidx = zeros(1,1,M); % index for each alternative
Cov_Noise_Inv_y_mu = zeros(1,N0,M); % for posterior mean and covariance 
Cov_Noise_Inv = cell(1,M); % for posterior mean and covariance
Cov_X2_V = cell(1,M); % for calculating regret in matrix form
x_indep_Cov_X_V = cell(1,M); % for calculating logIKG
for i = 1:M
    Cov_Noise_Inv{i} = zeros(N0,N0);
    Cov_X2_V{i} = zeros(J2,N0);
    x_indep_Cov_X_V{i} = zeros(J,N0);
end
N0cell = N0*ones(1,M);
x_indep_mu_prior = mu_prior(-M,X);
x_indep_mu = x_indep_mu_prior;

N1 = 2*C;

% start IKG
tStart = tic;

n = 1;
v0n = v0_pool(:,:,n);
% decide the sample location using IKG
[a,v,logIKGmax] = SGA_ave_bat_UV(IKG_pra,V,Vidx,Cov_Noise_Inv,Cov_Noise_Inv_y_mu,X,v0n,x_indep_mu,[],varPra,var_X,var_Cov_Inv_y_M);
consumption = cost(a,v,cost_type);

while consumption <= C
    % take sample at alternative a, location v
    Vidx(1,1,a) = Vidx(1,1,a) + 1;
    idx = Vidx(1,1,a);
    V(:,idx,a) = v;
    noise_var_a_v = noise_var(a,v); % simulation noise variance at this point
    yn(1,idx,a) = griewank_rev(a,v) + randn * sqrt(noise_var_a_v); % sample at this point; normal error
    mun(1,idx,a) = mu_prior(a,v); % prior mean at this point;

    if idx > 1
        %%% assume no covariance for noise
        cv = cov_vector (v, V(:,1:idx-1,a), theta(:,1,a), tau2(:,1,a))';
        A_inv = Cov_Noise_Inv{a}(1:idx-1,1:idx-1);
        alpha_inv = 1 / (tau2(:,1,a) + noise_var_a_v - cv'*A_inv*cv); 
        Cov_Noise_Inv{a}(1:idx,1:idx) = ...
            [A_inv + A_inv*(cv*cv')*A_inv*alpha_inv, -A_inv*cv*alpha_inv; ...
             -cv'*A_inv*alpha_inv,                   alpha_inv         ];
    else
        Cov_Noise_Inv{a}(1,1) = 1 / (tau2(:,1,a) + noise_var_a_v);
    end
    Cov_Noise_Inv_y_mu(1,1:idx,a) = (yn(1,1:idx,a) - mun(1,1:idx,a)) * Cov_Noise_Inv{a}(1:idx,1:idx); % store as row vector, more convenient    
    
    % compute regret
    Cov_X2_V{a}(:,idx) = cov_vector(v, X2, theta(:,1,a), tau2(:,1,a))';
    mu_i_n_x(a,:) = mu_i_n_x0(a,:) + Cov_Noise_Inv_y_mu(1,1:idx,a) * Cov_X2_V{a}(:,1:idx)'; 
    
    [~,chosen_index] = max(mu_i_n_x);
    index = sub2ind(size(theta_true),chosen_index,1:J2);
    true_chosen = theta_true(index); 
    regret_cost(:,n+1) = [mean(true_max - true_chosen); consumption];
    
    % update x_indep
    x_indep_Cov_X_V{a}(:,idx) = cov_vector(v, X, theta(:,1,a), tau2(:,1,a))';
    x_indep_mu(:,a) = x_indep_mu_prior(:,a) +  x_indep_Cov_X_V{a}(:,1:idx) * Cov_Noise_Inv_y_mu(1,1:idx,a)'; 
    
    % increase space
    if idx >= N0
        N0 = N0 + 500;

        store_temp = V;
        V = zeros(d,N0,M);
        V(:,1:idx,:) = store_temp;

        store_temp = yn;
        yn = zeros(1,N0,M);
        yn(1,1:idx,:) = store_temp;
        
        store_temp = mun;
        mun = zeros(1,N0,M);
        mun(1,1:idx,:) = store_temp;
        
        store_temp = Cov_Noise_Inv_y_mu;
        Cov_Noise_Inv_y_mu = zeros(1,N0,M);  
        Cov_Noise_Inv_y_mu(1,1:idx,:) = store_temp;
    end
    
    if idx >= N0cell(a)
        % increase space only for alternative a
        N0cell(a) = N0cell(a) + 500;
        
        store_temp = Cov_Noise_Inv{a};
        Cov_Noise_Inv{a} = zeros(N0cell(a),N0cell(a));
        Cov_Noise_Inv{a}(1:idx,1:idx) = store_temp;
        
        store_temp = Cov_X2_V{a};
        Cov_X2_V{a} = zeros(J2,N0cell(a));
        Cov_X2_V{a}(:,1:idx) = store_temp;
        
        store_temp = x_indep_Cov_X_V{a};
        x_indep_Cov_X_V{a} = zeros(J,N0cell(a));
        x_indep_Cov_X_V{a}(:,1:idx) = store_temp;
    end
    
    if n >= N1
        N1 = N1 + 500;
        v0_pool = cat(3,v0_pool, rand(d,M,500)*10);
        regret_cost = [regret_cost zeros(2,500)];
    end
    
    if mod(n,100) == 0
        fprintf('n=%d, consumption=%6.2f out of %6.2f, time= %6.2f seconds \n', n, consumption, C, toc(tStart));
    end
    
    % next round
    n = n + 1;
    v0n = v0_pool(:,:,n);
    % decide the sample location using IKG
    [a,v,logIKGmax] = SGA_ave_bat_UV(IKG_pra,V,Vidx,Cov_Noise_Inv,Cov_Noise_Inv_y_mu,X,v0n,x_indep_mu,x_indep_Cov_X_V,varPra,var_X,var_Cov_Inv_y_M);
    consumption = consumption + cost(a,v,cost_type);    
end
regret_cost_L{ell} = regret_cost(:,1:n);
end
delete(gcp)

% generate regret plot for L replications
regret_plot = zeros(L,C+1);
for ell = 1:L
    for c = 0:C
        idx = find(regret_cost_L{ell}(2,:) > c,1);
        if ~isempty(idx)
            regret_plot(ell,c+1) = regret_cost_L{ell}(1,idx-1);
        else
            regret_plot(ell,c+1) = regret_cost_L{ell}(1,end);
        end
    end
end

% save data
save(['IKGwSGA_UV_' num2str(d) 'D_density' num2str(density_type) '_cost' num2str(cost_type) '.mat'],'regret_cost_L','regret_plot');