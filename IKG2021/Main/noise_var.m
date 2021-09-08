function [lambada, lambada_d] = noise_var(i,x)

% constant sampling variance
lambada = 0.01 * ones(1,size(x,2));
% lambada = 0 * ones(1,size(x,2));   % for debug

if nargout > 1
    lambada_d = zeros(size(x));
end

end