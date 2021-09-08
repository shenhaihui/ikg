function X = rand_c(d,J,density_type)

if density_type == 1 % uniform on [0,10]^d
    X = rand(d,J)*10;
else % density_type == 2, independent N(0,sd^2), truncated on [0,10]^d
    sd = 4;
    X1 = trandn(zeros(d*J,1),10*ones(d*J,1)/sd); % truncated from standard normal
    X = reshape(sd*X1,d,[]); 
end

end

