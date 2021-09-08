function [noise_var_i_x, noise_var_d] = noise_var_est(i,x,varPra,var_X,var_Cov_Inv_y_M)

Cov_x_X = cov_vector(x, var_X, varPra.var_theta, varPra.var_tau2);
noise_var_i_x = Cov_x_X * var_Cov_Inv_y_M(:,i);

if nargout > 1
    noise_var_d = - diag(varPra.var_theta) * (repmat(x,1,varPra.k) - var_X) * diag(Cov_x_X) * var_Cov_Inv_y_M(:,i);
end

end