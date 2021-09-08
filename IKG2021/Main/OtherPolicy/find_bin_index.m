function loc_x = find_bin_index(x,bin_M)

x_grid = linspace(0,10,bin_M+1);
x_grid(1) = -10^(-6);

loc_x = find(x_grid >= x(1),1)-1;
if length(x) > 1
    for i = 2:length(x)
        loc_x = loc_x + (find(x_grid >= x(i),1)-1-1) * bin_M^(i-1);
    end
end
end