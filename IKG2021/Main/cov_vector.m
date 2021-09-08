function cv = cov_vector(x, X, theta, tau2)
% x - new point (d x 1), d dimension
% X - design points (d x k), d dimension, k points
% theta - factor of Gaussian covariance (d x 1 x M)
% tau2 - variance of radom field

n = size(X,2);
M = size(tau2,3); % M - alternative number

n1 = size(x,2);
if n1 == 1
    % Gaussian Covariance
    cv = repmat(tau2,1,n) .* exp(-sum(repmat(theta,1,n) .* (repmat(x,1,n,M) - X).^2,1));
else
    x2 = permute(x,[1 4 3 2]);
    cv = repmat(tau2,1,n,1,n1) .* exp(-sum(repmat(theta,1,n,1,n1) .* (repmat(x2,1,n,M,1) - repmat(X,1,1,1,n1)).^2,1));
end
    
end