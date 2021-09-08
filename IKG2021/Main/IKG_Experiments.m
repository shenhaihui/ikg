clear,clc; close all
addpath ./OtherPolicy

%% P1 & P2: 
% Constant sampling cost; density function (1) & (2); comparison with other
% three policies
cost_type = 1;
density_type = 1; % 1 or 2

for d = 1:2:7 % dimension of x
    SGA_m = 20*d;       % batch size
    SGA_K = 100*d;      % iteration number
    SGA_b1 = 200*d;     % b1  
    SGA_b2 = 0.7;       % b2; in SGA, b_k = b1/(k^(b2))
    parpool_n = 6;  % number of parallel works in computation
    
    if d == 1
        C = 50;
        if density_type == 1
            bin_M = 6;
        else
            bin_M = 8;
        end
    elseif d == 3
        C = 500;
        bin_M = 5;
    elseif d == 5
        C = 1000;
        if density_type == 1
            bin_M = 3;
        else
            bin_M = 2;
        end        
    elseif d == 7
        C = 1500;
        parpool_n = 4;
        bin_M = 1;
    end
    
    % IKG
    regret_plot = IKGwSGA(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2);
    figure, plot(0:C, mean(regret_plot), 'g-'); hold on;
    % IKGwRC
    regret_plot = IKGwRC(density_type,cost_type,parpool_n,d,C);
    plot(0:C, mean(regret_plot), 'b--');
    % BSE
    regret_plot = BSE(density_type,cost_type,parpool_n,d,C,bin_M);
    plot(0:C, mean(regret_plot), 'm:');
    % PRS
    regret_plot = PRS(density_type,cost_type,parpool_n,d,C);
    plot(0:C, mean(regret_plot), 'r-.');
    legend('IKG','IKGwRC','BSE','PRS');
    title(['d=', num2str(d)]);
end


%% P3: 
% Varying sampling cost; density function (1)
% IKG
% IKG-CI: ignore variations in the sampling cost at different locations
% and mistakenly uses the unit sampling cost when implementing the IKG,
% but the actual sampling consumption follows c_i(x)
cost_type = 2;
density_type = 1;

for d = 1:2:7 % dimension of x
    SGA_m = 20*d;       % batch size
    SGA_K = 100*d;      % iteration number
    SGA_b1 = 200*d;     % b1  
    SGA_b2 = 0.7;       % b2; in SGA, b_k = b1/(k^(b2))
    parpool_n = 6;  % number of parallel works in computation
    
    if d == 1
        C = 50;
    elseif d == 3
        C = 500;
    elseif d == 5
        C = 1000;   
    elseif d == 7
        C = 1500;
        parpool_n = 4;
    end
    
    % IKG
    regret_plot = IKGwSGA(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2);
    figure, plot(0:C, mean(regret_plot), 'g-'); hold on;
    E = 2.7564 * sqrt(var(regret_plot)/size(regret_plot,1)); % 99% CI
    plot(0:C, mean(regret_plot) + E, 'g:');
    h1 = plot(0:C, mean(regret_plot) - E, 'g:');
    % IKG-CI
    regret_plot = IKGwSGA_CI(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2);
    plot(0:C, mean(regret_plot), 'r--');
    E = 2.7564 * sqrt(var(regret_plot)/size(regret_plot,1)); % 99% CI
    plot(0:C, mean(regret_plot) + E, 'r:');
    plot(0:C, mean(regret_plot) - E, 'r:');
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('IKG','+- 99% CI','IKG-CI','+- 99% CI');
    title(['d=', num2str(d)]);
end