function g_x = h_c_partial_x_X(V,Vidx,theta,tau2,Cov_Noise_Inv,X,i,x,x_indep_mu,x_indep_Cov_X_V,cost_type)

% v - a random sample from domain
% i - considered alternative
% x - derivative of h_i_n(v,x) with respect to x

%%% compute mu and cov_X_V out of this function

% M = size(V,3); % M - alternative number
% [d,~,M] = size(V); % d - dimension, M - alternative number
[d,J] = size(X);

[~,noise_var_d] = noise_var(i,x); 


idx = Vidx(1,1,i);

if idx == 0
    k_i_n_x_x = cov_vector(x, x, theta(:,1,i), tau2(:,1,i));
    k_i_n_x_X = cov_vector(x, X, theta(:,1,i), tau2(:,1,i));
    
    lam = k_i_n_x_x + noise_var(i,x);
    sigma_tilde_X = k_i_n_x_X * lam^(-1/2);    
    
    partial_sigma = 2 * lam^(-1/2) * (repmat(cov_vector(x, X, theta(:,1,i), tau2(:,1,i)),d,1) ...
        .* repmat(theta(:,1,i),1,J) .* (X - repmat(x,1,J)))...
        -1/2 * lam^(-3/2) * repmat(k_i_n_x_X,d,1) .* repmat(noise_var_d,1,J);       
    
else
    cv = cov_vector(x, V(:,1:idx,i), theta(:,1,i), tau2(:,1,i));
    k_i_n_x_x = cov_vector(x, x, theta(:,1,i), tau2(:,1,i)) - ...
        cv * Cov_Noise_Inv{i}(1:idx,1:idx) * cv';
    
    if k_i_n_x_x < 0 % may be caused by numerical error!
        k_i_n_x_x = -k_i_n_x_x;
    end
    
    k_i_n_x_X = cov_vector(x, X, theta(:,1,i), tau2(:,1,i)) - ...
        (x_indep_Cov_X_V{i}(:,1:idx) * Cov_Noise_Inv{i}(1:idx,1:idx) * cov_vector(x, V(:,1:idx,i), theta(:,1,i), tau2(:,1,i))')';      
    
    lam = k_i_n_x_x + noise_var(i,x);
    sigma_tilde_X = k_i_n_x_X * lam^(-1/2);        
    
    A = diag(theta(:,1,i)) * (repmat(x,1,idx) - V(:,1:idx,i)) ...
          * diag(cov_vector(x, V(:,1:idx,i), theta(:,1,i), tau2(:,1,i))) ...
          * Cov_Noise_Inv{i}(1:idx,1:idx);            
      
    partial_sigma = 2 * lam^(-1/2) * (repmat(cov_vector(x, X, theta(:,1,i), tau2(:,1,i)),d,1) ...
        .* repmat(theta(:,1,i),1,J) .* (X - repmat(x,1,J)) + A * x_indep_Cov_X_V{i}(:,1:idx)')...
        -1/2 * lam^(-3/2) * repmat(k_i_n_x_X,d,1) .* repmat((4 * A * cov_vector(x, V(:,1:idx,i), theta(:,1,i), tau2(:,1,i))' + noise_var_d),1,J);         
end


I = abs(sigma_tilde_X) < 1e-10;
sigma_tilde_X(I) = 1;

mu_i = x_indep_mu(:,i);
mu_no_i = x_indep_mu;
mu_no_i(:,i)=[];
beta = abs(max(mu_no_i,[],2) - mu_i)' ./ sigma_tilde_X;
h_p_x = repmat(sign(sigma_tilde_X) .* normpdf(beta),d,1) .* partial_sigma;
h_p_x(:,I) = 0;

% average of h_p_x
h_p_x = mean(h_p_x,2);

% consider the cost
[c, d_c_x] = cost(i,x,cost_type);
if d_c_x == 0
    g_x = h_p_x / c;
else
    neg_logIKG = neg_logIKG_i_x(V,Vidx,theta,tau2,Cov_Noise_Inv,X,i,x,x_indep_mu,x_indep_Cov_X_V,cost_type);
    h2c = exp(-neg_logIKG);
    g_x = (h_p_x - h2c * d_c_x) / c;
end

end