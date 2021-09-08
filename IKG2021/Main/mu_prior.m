function mu0 = mu_prior(i,x)

% zero prior mean

if i > 0
    mu0 = 0;
else % if i == -M, output for all i=1,...M and all points in x
    mu0 = zeros(size(x,2),-i);
end

end
