function [a,v,logIKGmax] = SAA(pra,V,Vidx,Cov_Noise_Inv,X,x0,x_indep_mu,x_indep_Cov_X_V,x_indep_mu_MB,x_indep_Cov_X_V_MB,v_MB)

lb = pra.lb; % lower bound
ub = pra.ub; % upper bound
theta = pra.theta;
tau2 = pra.tau2;

cost_type = pra.cost_type;

[d,~,M] = size(V); % d - dimension; M - alternative number

% options = optimoptions('fmincon','Display','iter','GradObj','on');
options = optimoptions('fmincon','Display','off','GradObj','on');
for i = 1:M
    x_initial = x0(:,i);
    x0(:,i) = fmincon(@(x) neg_logIKG_i_x(V,Vidx,theta,tau2,Cov_Noise_Inv,v_MB,i,x,x_indep_mu_MB,x_indep_Cov_X_V_MB,cost_type),...
        x_initial,[],[],[],[],lb,ub,[],options);         
end

neg_logIKG = zeros(1,M);
for i = 1:M        
    neg_logIKG(i) = neg_logIKG_i_x(V,Vidx,theta,tau2,Cov_Noise_Inv,X,i,x0(:,i),x_indep_mu,x_indep_Cov_X_V,cost_type);
end

[logIKGmax,a] = max(-neg_logIKG);
v = x0(:,a);
end