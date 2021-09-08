function [neg_logIKG, neg_d_logIKG] = neg_logIKG_i_x_UV(V,Vidx,theta,tau2,Cov_Noise_Inv,X,i,x,x_indep_mu,x_indep_Cov_X_V,cost_type,varPra,var_X,var_Cov_Inv_y_M)

% for given i (alternative) and x (point), compute neg_logIKG
% based on SSA with sample X

%%% compute mu and cov_X_V out of this function

J = size(X,2); % number of sampled points to estimate integral
idx = Vidx(1,1,i);
if idx == 0
    k_i_n_x_x = cov_vector(x, x, theta(:,1,i), tau2(:,1,i));
    lam = sqrt(k_i_n_x_x + noise_var_est(i,x,varPra,var_X,var_Cov_Inv_y_M));
    sigma_tilde_abs = abs(cov_vector(x, X, theta(:,1,i), tau2(:,1,i)))' / lam;
else
    cv = cov_vector(x, V(:,1:idx,i), theta(:,1,i), tau2(:,1,i));
    k_i_n_x_x = cov_vector(x, x, theta(:,1,i), tau2(:,1,i)) - ...
        cv * Cov_Noise_Inv{i}(1:idx,1:idx) * cv';
    
    if k_i_n_x_x < 0 % caused by numerical error!
        k_i_n_x_x = -k_i_n_x_x;
    end
    
    lam = sqrt(k_i_n_x_x + noise_var_est(i,x,varPra,var_X,var_Cov_Inv_y_M));
    
    sigma_tilde_abs = cov_vector(x, X, theta(:,1,i), tau2(:,1,i))' - ...
        x_indep_Cov_X_V{i}(:,1:idx) * Cov_Noise_Inv{i}(1:idx,1:idx) * cov_vector(x, V(:,1:idx,i), theta(:,1,i), tau2(:,1,i))';
    sigma_tilde_abs = abs(sigma_tilde_abs) / lam;
end


% --- compute g ---
mu_i = x_indep_mu(:,i);
mu_no_i = x_indep_mu;
mu_no_i(:,i)=[];
Delta_abs = abs(max(mu_no_i,[],2) - mu_i);
sigma_tilde_abs_0 = (sigma_tilde_abs < 1e-10);
sigma_tilde_abs(sigma_tilde_abs_0) = [];
Delta_abs(sigma_tilde_abs_0) = [];
u = Delta_abs ./ sigma_tilde_abs;
p = 1/J * ones(size(sigma_tilde_abs));
g_temp = log(1/sqrt(2*pi)) +  log(p .* sigma_tilde_abs) - u.^2/2;     
I = u > 20;
u_r_appro = u .* u ./ (u.^2+1);
u(I) = 0;
u_r_exact = u .* normcdf(-u) ./ normpdf(u); % for large theta, it may give NaN
u_r = u_r_appro .* I + u_r_exact .* (~I);
g = g_temp + log1p(-u_r); % log1p() may give -Inf

% --- compute logIKG ---
if isempty(g)
    logIKG = -inf;
else
    gmax = max(g);
    if gmax > -inf
        logIKG = gmax + log(sum(exp(g-gmax)));
    else
        logIKG = -inf;
    end
end

% neg_logIKG = - logIKG;

% consider the cost
neg_logIKG = - logIKG + log(cost(i,x,cost_type));



%%% supply gradient
if nargout > 1 
    g_x = h_c_partial_x_X(V,Vidx,theta,tau2,Cov_Noise_Inv,X,i,x,x_indep_mu,x_indep_Cov_X_V,cost_type);
    
    if logIKG > -inf
        neg_d_logIKG = - g_x / exp(-neg_logIKG); 
    else
        neg_d_logIKG = - g_x * 1000;
        fprintf('logIKG == -inf\n')
    end
end

end