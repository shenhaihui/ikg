clear,clc; close all
addpath ../../Main

% P1: Constant sampling cost; density function (1)
cost_type = 1;
density_type = 1;
C = 100;

%% Tune the sample size for IKGwSAA
% change the value of d within {1,2,...,7}
d = 1;  % dimension of x
SGA_m = 20*d;       % batch size
SGA_K = 100*d;      % iteration number
SGA_b1 = 200*d;     % b1  
SGA_b2 = 0.7;       % b2; in SGA, b_k = b1/(k^(b2))
parpool_n = [];  % do not use parallel computation, to better compare time

% IKGwSGA
regret_plot = IKGwSGA(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2);
figure, plot(0:C, mean(regret_plot), 'g-'); hold on;

% IKGwSAA
% For d = 1,2,...,7, try different value of sample_size_ratio so that the
% resulting OC curve from IKGwSAA is roughly the same to that from IKGwSGA.
sample_size_ratio = 1/16; 
regret_plot = IKGwSAA(density_type,cost_type,parpool_n,d,C,sample_size_ratio);
plot(0:C, mean(regret_plot), 'r-');


%% Run the entire loop for d = 1,2,...,7
% After trial, the found values of sample_size_ratio for different d are as
% follows: d==1~3, 1/16; d==4, 1/4; d==5, 1; d==6, 2; d==7, 4
for d = 1:7 % dimension of x
    SGA_m = 20*d;       % batch size
    SGA_K = 100*d;      % iteration number
    SGA_b1 = 200*d;     % b1  
    SGA_b2 = 0.7;       % b2; in SGA, b_k = b1/(k^(b2))
    parpool_n = [];  % do not use parallel computation, to better compare time

    % IKGwSGA
    regret_plot = IKGwSGA(density_type,cost_type,parpool_n,d,C,SGA_m,SGA_K,SGA_b1,SGA_b2);
    figure, plot(0:C, mean(regret_plot), 'g-'); hold on;

    % IKGwSAA
    if d <= 3
        sample_size_ratio = 1/16; 
    elseif d == 4
        sample_size_ratio = 1/4;
    elseif d == 5
        sample_size_ratio = 1;
    elseif d == 6
        sample_size_ratio = 2;
    elseif d == 7
        sample_size_ratio = 4;
    end
    regret_plot = IKGwSAA(density_type,cost_type,parpool_n,d,C,sample_size_ratio);
    plot(0:C, mean(regret_plot), 'r-');
end

%% Plot result
d = 1:7;
N1 = 2000*d.^2;
N2 = N1 .* [1/16, 1/16, 1/16, 1/4, 1, 2, 4];
figure, plot(d, N1, 'g-o'); hold on;
plot(d, N2, 'r-s');
xlabel('d');
ylabel('Sample Size');
legend('N_1(d)','N_2(d)');

T1 = zeros(1,7);
T2 = zeros(1,7);
for d = 1:7
    load(['IKGwSGA_' num2str(d) 'D_density1_cost1.mat']);
    T1(d) = mean(time_L);
    load(['IKGwSAA_' num2str(d) 'D_density1_cost1.mat']);
    T2(d) = mean(time_L);
end
d = 1:7;
figure, plot(d, T1, 'g-o'); hold on;
plot(d, T2, 'r-s');
xlabel('d');
ylabel('Computation Time (sec.)');
legend('T_1(d)','T_2(d)');