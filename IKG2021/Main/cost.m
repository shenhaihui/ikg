function [c, d_c_x] = cost(i,x,cost_type)

if cost_type == 1 % cost is identical;
    c = 1;
    if nargout > 1
        d_c_x = 0;
    end    
else 
    d = size(x,1);
    c = 2^(3-i) * (1 + sum((x-5).^2) / (10*d));
    if nargout > 1
        d_c_x = 2^(3-i) * 2 * (x-5) / (10*d);
    end  
end

end