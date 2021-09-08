function [lambada, lambada_d] = noise_var(i,x)

%%% constant sampling variance, 0.01
lambada = 0.01 * ones(1,size(x,2));
% lambada = 0 * ones(1,size(x,2));   % for debug

if nargout > 1
    lambada_d = zeros(size(x));
end

% %%% varying  sampling variance
% [d, m]=size(x); % d - dimension, m - number of points
% lambada = 0.01 * (1.5^(d-1) + griewank_rev(i,x));
% 
% if nargout > 1
%     cos_temp = cos(x ./ repmat(sqrt(i*(1:d)'),1,m));
%     cos_leave_k = zeros(d,m);
%     for k = 1:d
%         cos_temp2 = cos_temp;
%         cos_temp2(k,:) = ones(1,m);
%         cos_leave_k(k,:) = prod(cos_temp2,1);
%     end
%     lambada_d = 0.01 * (x/2000 + 1.5^(d-1) * ...
%         sin(x ./ repmat(sqrt(i*(1:d)'),1,m)) ./ repmat(sqrt(i*(1:d)'),1,m) .* cos_leave_k);
% end

end