clear,clc; %close all
addpath ../../Main

% P1
% Constant sampling cost; density function (1)
cost_type = 1;
density_type = 1;

for d = 1:2:3  % dimension of x
    SGA_m = 20*d;       % batch size
    SGA_K = 100*d;      % iteration number
    SGA_b1 = 200*d;     % b1  
    SGA_b2 = 0.7;       % b2; in SGA, b_k = b1/(k^(b2))
    parpool_n = 6;  % number of parallel works in computation

    if d == 1
        C = 50;
    elseif d == 3
        C = 500;
    end

    % with estimated variance surface
    regret_plot = IKGwSGA_UV(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2);
    figure, plot(0:C, mean(regret_plot), 'b--'); hold on;

    % with true variance surface
    regret_plot = IKGwSGA(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2);
    plot(0:C, mean(regret_plot), 'g-');
    legend('Estimated var.','True var.');
    title(['d=', num2str(d)]);
end