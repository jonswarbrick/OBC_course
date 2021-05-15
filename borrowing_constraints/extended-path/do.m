%% Solves SOE model with a borrowing constraint using dynare extended-path method
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

%% With OBC
dynare soe_borrowing_constraint
clear;

%% Without OBC
dynare soe
clear;

%% Cleanup
dynare_cleanup(  );

%% Plot results and report moments
burn_in = 1000;
plot_start = 5000;
plot_end = 5500;

load('results/simulated_time_series.mat')
sim1 = sims; offset1 = offset;
load('results/simulated_time_series_unbounded.mat')
sim2 = sims; offset2 = offset;
red_1 = [153 0 51] / 255;

h = figure;
subplot(2,2,1); plot(log(sim1.c(plot_start:plot_end)./offset1.c),'k-','LineWidth',1); title('log(c)'); axis tight; hold on;
                plot(log(sim2.c(plot_start:plot_end)./offset2.c),'Color',red_1,'LineWidth',1); 
subplot(2,2,2); plot(log(sim1.h(plot_start:plot_end)./offset1.h),'k-','LineWidth',1); title('log(h)'); axis tight; hold on;
                plot(log(sim2.h(plot_start:plot_end)./offset2.h),'Color',red_1,'LineWidth',1);  
subplot(2,2,3); plot(sim1.b(plot_start:plot_end)./sim1.c(plot_start:plot_end),'k-','LineWidth',1); title('b / c'); axis tight; hold on;
                plot(sim2.b(plot_start:plot_end)./sim2.c(plot_start:plot_end),'Color',red_1,'LineWidth',1);
subplot(2,2,4); plot(sim1.mu(plot_start:plot_end)./sim1.c(plot_start:plot_end),'k-','LineWidth',1); title('\mu / c'); axis tight; hold on;
                plot(sim2.mu(plot_start:plot_end)./sim2.c(plot_start:plot_end),'Color',red_1,'LineWidth',1); 
h.Position = [ 100 , 100 , 800 , 400 ];
print(h,'ep_sims','-depsc')
legend('Borrowing Constraint','No constraint')

% Metrics
moments(1,1) = mean(sim1.c(burn_in+1:end));
moments(2,1) = mean(sim1.h(burn_in+1:end));
moments(3,1) = mean(sim1.b(burn_in+1:end)/mean(sim1.c(burn_in+1:end)));
moments(1,2) = std(sim1.c(burn_in+1:end)/mean(sim1.c(burn_in+1:end)));
moments(2,2) = std(sim1.h(burn_in+1:end)/mean(sim1.h(burn_in+1:end)));
moments(3,2) = std(sim1.b(burn_in+1:end)/mean(sim1.c(burn_in+1:end)));
moments(1,3) = skewness(sim1.c(burn_in+1:end));
moments(2,3) = skewness(sim1.h(burn_in+1:end));
moments(3,3) = skewness(sim1.b(burn_in+1:end));

disp('**-- Borrowing constraints model --** ')
disp( table( moments(:,1) , moments(:,2) , moments(:,3) , ...
          'VariableNames',{'Mean','StandardDeviation','Skewness'},...
          'RowNames',{'Consumption';'Hours';'Bonds'}) )


con_binding = zeros(size(sim1.b(burn_in+1:end)));
con_binding(sim1.b(burn_in+1:end)+0.01<1e-6) = 1;
disp(['Constraint binds in ', num2str(100*mean(con_binding)),'% of periods'])

moments(1,1) = mean(sim2.c(burn_in+1:end));
moments(2,1) = mean(sim2.h(burn_in+1:end));
moments(3,1) = mean(sim2.b(burn_in+1:end)/mean(sim2.c(burn_in+1:end)));
moments(1,2) = std(sim2.c(burn_in+1:end)/mean(sim2.c(burn_in+1:end)));
moments(2,2) = std(sim2.h(burn_in+1:end)/mean(sim2.h(burn_in+1:end)));
moments(3,2) = std(sim2.b(burn_in+1:end)/mean(sim2.c(burn_in+1:end)));
moments(1,3) = skewness(sim2.c(burn_in+1:end));
moments(2,3) = skewness(sim2.h(burn_in+1:end));
moments(3,3) = skewness(sim2.b(burn_in+1:end));

disp('**-- Unbounded --** ')
disp( table( moments(:,1) , moments(:,2) , moments(:,3) , ...
          'VariableNames',{'Mean','StandardDeviation','Skewness'},...
          'RowNames',{'Consumption';'Hours';'Bonds'}) )
