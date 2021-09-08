function [a,v,logIKGmax] = NoSGA(pra,V,Vidx,Cov_Noise_Inv,X,x0,x_indep_mu,x_indep_Cov_X_V)

theta = pra.theta;
tau2 = pra.tau2;

cost_type = pra.cost_type;

[d,~,M] = size(V); % d - dimension; M - alternative number

neg_logIKG = zeros(1,M);
for i = 1:M        
    neg_logIKG(i) = neg_logIKG_i_x(V,Vidx,theta,tau2,Cov_Noise_Inv,X,i,x0(:,i),x_indep_mu,x_indep_Cov_X_V,cost_type);
end

[logIKGmax,a] = max(-neg_logIKG);
v = x0(:,a);
end

