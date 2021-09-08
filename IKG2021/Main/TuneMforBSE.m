clear,clc
addpath ./OtherPolicy

% cost_type = 1;
% bin_M_range = [1 2 3 4 5 6 7 8 9 10];
% BSE_regret_reduction = zeros(10,8);
% 
% for i = 1:10
%     bin_M = bin_M_range(i);
%     k = 0;
%     for density_type = 1:2
%         parpool_n = 6;
%         for d = 1:2:7
%             if d == 1
%                 C = 50;
%             elseif d == 3
%                 C = 500;
%             elseif d == 5
%                 C = 1000;
%             elseif d == 7
%                 C = 1500;
%                 parpool_n = 4;
%             end
%             k = k+1;
%             regret_plot = BSE(density_type,cost_type,parpool_n,d,C,bin_M);
%             regret_ave = mean(regret_plot);
%             BSE_regret_reduction(i,k) = regret_ave(end) / regret_ave(1);
%         end
%     end
% end
% save('BSE_regret_reduction.mat','BSE_regret_reduction');

load('BSE_regret_reduction.mat')
k = 0;
for density_type = 1:2
    for d = 1:2:7
        k = k+1;
        [~,M_best] = min(BSE_regret_reduction(:,k));
        fprintf('density_type=%d, dim=%d, M_best=%d \n', density_type, d, M_best);
    end
end