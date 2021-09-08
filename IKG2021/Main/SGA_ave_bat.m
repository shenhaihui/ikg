% -- input -- 
% Xn - visited points, for each slice (alternative), each colum is a point
% XYidx - index for each alternative
% theta, tau2 -  parameters for Gaussian Covariance
% Cov_Noise_Inv, Cov_Noise_Inv_y_mu - for posterior mean and covariance
% X - sample points, each colum is a point
% lb, ub - bound while searching best x

% b - control step size
% SGD_N - iteration number

% -- output --
% (a,v,loghmax) - alternative, location, corresponding logh

function [a,v,logIKGmax] = SGA_ave_bat(pra,V,Vidx,Cov_Noise_Inv,Cov_Noise_Inv_y_mu,X,x0,x_indep_mu,x_indep_Cov_X_V)

lb = pra.lb; % lower bound
ub = pra.ub; % upper bound
b1 = pra.SGA_b1; % control step size, b1/k^b2
b2 = pra.SGA_b2;
K = pra.SGA_K; % iteration number
m = pra.SGA_m; % mini batch size
theta = pra.theta;
tau2 = pra.tau2;

density_type = pra.density_type;
cost_type = pra.cost_type;

[d,~,M] = size(V); % d - dimension; M - alternative number
LB = repmat(lb,1,M);
UB = repmat(ub,1,M);
idx_max = max(Vidx);

x0_temp = zeros(d,M,K);
for k = 1:K
    v_MB = rand_c(d,m,density_type); % Mini Batch 
    [x_indep_mu_MB, x_indep_Cov_X_V_MB] = x_indep_comp(V,idx_max,Cov_Noise_Inv_y_mu,v_MB,M,theta,tau2);
    for i = 1:M
        g_x = h_c_partial_x_X(V,Vidx,theta,tau2,Cov_Noise_Inv,v_MB,i,x0(:,i),x_indep_mu_MB,x_indep_Cov_X_V_MB,cost_type); % averaged gradient over v_MB on x0(:,i)
        x0(:,i) = x0(:,i) + b1/k^b2 * g_x;
    end 
    % project x back into the domain of x; only for rectangle domain
    I1 = x0 < LB;
    I2 = x0 > UB;
    x0 = LB.*I1 + UB.*I2 + x0.*(~I1.*~I2);   

    x0_temp(:,:,k) = x0;          
end
x0_stop = mean(x0_temp(:,:,3/4*K:K),3);

neg_logIKG = zeros(1,M);
for i = 1:M        
    neg_logIKG(i) = neg_logIKG_i_x(V,Vidx,theta,tau2,Cov_Noise_Inv,X,i,x0_stop(:,i),x_indep_mu,x_indep_Cov_X_V,cost_type);
end

[logIKGmax,a] = max(-neg_logIKG);
v = x0_stop(:,a);
end


function [mu, Cov_X_V] = x_indep_comp(V,idx_max,Cov_Noise_Inv_y_mu,X,M,theta,tau2)
J = size(X,2);
if idx_max == 0
    mu = mu_prior(-M,X);
    Cov_X_V = [];
else
    Cov_X_V = cov_vector(X,V(:,1:idx_max,:),theta,tau2);
    Cov_X_V = permute(Cov_X_V,[4,2,3,1]);     
    % if out of memory, use the following instead
%     Cov_X_V = zeros(J,idx_max,M);
%     for j = 1:J
%         xj = X(:,j);
%         Cov_X_V(j,:,:) = cov_vector(xj,V(:,1:idx_max,:),theta,tau2);
%     end    
    
    mu_temp = sum(Cov_X_V .* repmat(Cov_Noise_Inv_y_mu(1,1:idx_max,:),J,1),2);
    mu = mu_prior(-M,X) + permute(mu_temp,[1 3 2]);
    
    % store in cell type
    Cov_X_V = num2cell(Cov_X_V,[1 2]);
end
end
