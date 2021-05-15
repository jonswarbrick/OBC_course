%% Solves SOE model using RISE
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
path('../../functions/RISE_toolbox-master',path);
path('../../functions',path);
rise_startup();

% Create the model structure from the model file
m=rise('soe_borrowing_constraint' );

% Solve the policy functions
[ms,retcode]=solve(m);
% m = set( m , 'solve_derivatives_type','automatic');
%[ms,retcode]=solve(m , 'solve_order' , 2 );

% Print out policy function
print_solution(ms);

%% Compute average IRFs
% girf=irf(ms,'irf_periods',100,'irf_shock_sign',-2,'irf_anticipate',false,...
%         'simul_honor_constraints_through_switch',true,'simul_honor_constraints',true,...
%         'irf_regime_specific',false,'irf_type','girf','irf_draws',500);
% 
% irfs.c.epsz = girf.epsz.c.data;
% irfs.b.epsz = girf.epsz.b.data;
% irfs.h.epsz = girf.epsz.h.data;
% save('results/rise_irfs','irfs');
% 
% figure;
% subplot( 1 , 3 , 1)
% plot( irfs.c.epsz(1:60),  'r-' , 'LineWidth', 2  );
% title('Consumption')
% grid on
% subplot( 1 , 3 , 2)
% plot( irfs.h.epsz(1:60),  'r-' , 'LineWidth', 2  );
% title('Hours')
% grid on
% subplot( 1 , 3 , 3)
% plot( irfs.b.epsz(1:60),  'r-' , 'LineWidth', 2  );
% title('Bonds')
% grid on
% set(gcf,'position',[200 200 580 130])
% set(findall(gcf,'-property','FontSize'),'FontSize',8)

%% Simulate time-series
H = 300000;
burn_in = 1000;
b_limit = -0.01;

rng(147);
[db,states,retcode] = simulate(m,'simul_periods',H,'simul_honor_constraints_through_switch',true, 'simul_honor_constraints',true , 'simul_burn' , 0 );

cons = db.c.data(burn_in+1:end);
hours = db.h.data(burn_in+1:end);
prod = db.z.data(burn_in+1:end);
bonds = db.b.data(burn_in+1:end);

% Metrics
moments(1,1) = mean(cons);
moments(2,1) = mean(hours);
moments(3,1) = mean(bonds/mean(cons));
moments(1,2) = std(cons/mean(cons));
moments(2,2) = std(hours/mean(hours));
moments(3,2) = std(bonds/mean(cons));
moments(1,3) = skewness(cons);
moments(2,3) = skewness(hours);
moments(3,3) = skewness(bonds);

format long

disp( table( moments(:,1) , moments(:,2) , moments(:,3) , ...
          'VariableNames',{'Mean','StandardDeviation','Skewness'},...
          'RowNames',{'Consumption';'Hours';'Bonds'}) )

con_binding = zeros(size(bonds));
con_binding(bonds-b_limit<1e-6) = 1;
disp(['Constraint binds in ', num2str(100*mean(con_binding)),'% of periods'])

format short
