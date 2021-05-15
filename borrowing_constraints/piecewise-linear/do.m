%% Solves SOE model using piecewise-linear method
% Uses Iacoviello and Guerrieri's 'OccBin' toolbox
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
path('../../functions/occbin',path);
path('../../functions',path);

% .mod filenames for two regimes
opts.modnam = 'soe_borrowing_constraint_slack';
opts.modnamstar = 'soe_borrowing_constraint_binding';

opts.constraint = 'b<b_limit';
opts.constraint_relax ='mu<0';

opts.simulate_girf = 0; % 1 or 0

%% IRF simulation

% Pick innovation for IRFs
irfshock =char('epsz');      % label for innovation for IRFs

shockscale = -2;
nperiods = 80;
shockssequence = zeros(nperiods,1);
shockssequence(2) = 1;
shockssequence = shockssequence * shockscale;
    
[zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
  solve_one_constraint(opts.modnam,opts.modnamstar,...
  opts.constraint, opts.constraint_relax,...
  shockssequence,irfshock,nperiods);

% Unpack simulations
for i=1:Mbase_.endo_nbr
  piecewise.irf_ss.(deblank(Mbase_.endo_names(i,:))).epsz = zdatass(i);
  if strcmp(deblank(Mbase_.endo_names(i,:)),'b')
      piecewise.irf_unbounded.(deblank(Mbase_.endo_names(i,:))).epsz = zdatalinear(2:61,i)' / zdatass(1);
      piecewise.irf.(deblank(Mbase_.endo_names(i,:))).epsz = zdatapiecewise(2:61,i)' / zdatass(1);
  else
      piecewise.irf_unbounded.(deblank(Mbase_.endo_names(i,:))).epsz = zdatalinear(2:61,i)' / zdatass(i);
      piecewise.irf.(deblank(Mbase_.endo_names(i,:))).epsz = zdatapiecewise(2:61,i)' / zdatass(i);
  end
end


save('results/irfs.mat','piecewise');

%% GIRF
if opts.simulate_girf
% Pick innovation for IRFs
opts.irfshock =char('epsz');
opts.replic = 200;
opts.irf_period = 80;
opts.burnin = 100;
opts.burnout = 20;
opts.shockscale = -2;

for mc=1:opts.replic
    shockscale = -2;
    nperiods = opts.burnin+opts.irf_period+opts.burnout;
    rng(mc);
    shockssequence = randn(nperiods,1);

    try
    [zdatalinear_1 zdatapiecewise_1 zdatass_1 oobase_1 Mbase_1  ] = ...
  solve_one_constraint(opts.modnam,opts.modnamstar,...
  opts.constraint, opts.constraint_relax,...
  shockssequence,opts.irfshock,nperiods);
    catch
        pause(2);
    [zdatalinear_1 zdatapiecewise_1 zdatass_1 oobase_1 Mbase_1  ] = ...
  solve_one_constraint(opts.modnam,opts.modnamstar,...
  opts.constraint, opts.constraint_relax,...
  shockssequence,opts.irfshock,nperiods);
    end
  
    shockssequence(opts.burnin+1) = shockssequence(opts.burnin+1)+shockscale;
    
    try
    [zdatalinear_2 zdatapiecewise_2 zdatass_2 oobase_2 Mbase_2  ] = ...
  solve_one_constraint(opts.modnam,opts.modnamstar,...
  opts.constraint, opts.constraint_relax,...
  shockssequence,opts.irfshock,nperiods);
    catch
        pause(2);
    [zdatalinear_2 zdatapiecewise_2 zdatass_2 oobase_2 Mbase_2  ] = ...
  solve_one_constraint(opts.modnam,opts.modnamstar,...
  opts.constraint, opts.constraint_relax,...
  shockssequence,opts.irfshock,nperiods);
    end
  
% Unpack simulations
for i=1:Mbase_1.endo_nbr
  piecewise.girf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz_mat(mc,:) = zdatalinear_2(opts.burnin+1:opts.burnin+opts.irf_period,i)'-zdatalinear_1(opts.burnin+1:opts.burnin+opts.irf_period,i)';
  piecewise.girf.(deblank(Mbase_1.endo_names(i,:))).epsz_mat(mc,:) = zdatapiecewise_2(opts.burnin+1:opts.burnin+opts.irf_period,i)'-zdatapiecewise_1(opts.burnin+1:opts.burnin+opts.irf_period,i)';
end

disp(mc)
end

% Unpack simulations
for i=1:Mbase_1.endo_nbr
  piecewise.girf_ss.(deblank(Mbase_1.endo_names(i,:))).epsz = zdatass(i);
  if strcmp(deblank(Mbase_1.endo_names(i,:)),'b')
      piecewise.girf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz = mean(piecewise.girf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(1);
      piecewise.girf.(deblank(Mbase_1.endo_names(i,:))).epsz = mean(piecewise.girf.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(1);
      piecewise.median_irf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz = median(piecewise.girf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(1);
      piecewise.median_girf.(deblank(Mbase_1.endo_names(i,:))).epsz = median(piecewise.girf.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(1);
  else
      piecewise.girf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz = mean(piecewise.girf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(i);
      piecewise.girf.(deblank(Mbase_1.endo_names(i,:))).epsz = mean(piecewise.girf.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(i);
      piecewise.median_irf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz = median(piecewise.girf_unbounded.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(i);
      piecewise.median_girf.(deblank(Mbase_1.endo_names(i,:))).epsz = median(piecewise.girf.(deblank(Mbase_1.endo_names(i,:))).epsz_mat) / zdatass(i);
  end
end

save('results/irfs.mat','piecewise');
end

%% Time-series simulation

% Pick innovation for IRFs
irfshock =char('epsz');      % label for innovation for IRFs

rng(147);
nperiods = 10000;
burn_in = 1000;
shockssequence = randn(nperiods,1);
    
[zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
  solve_one_constraint(opts.modnam,opts.modnamstar,...
  opts.constraint, opts.constraint_relax,...
  shockssequence,irfshock,nperiods);

% Unpack simulations
for i=1:Mbase_.endo_nbr
  piecewise_sim.sim_unbounded.(deblank(Mbase_.endo_names(i,:))) = zdatalinear(:,i)'+ zdatass(i)';
  piecewise_sim.sim.(deblank(Mbase_.endo_names(i,:))) = zdatapiecewise(:,i)' + zdatass(i)';
end

save('results/sims.mat','piecewise');

%% Figures and Moments
figure;

subplot(1,3,1)
plot(piecewise.irf_unbounded.c.epsz,'k-','LineWidth',2); hold on;
plot(piecewise.irf.c.epsz,'r-','LineWidth',2); 
title('Consumption')
subplot(1,3,2)
plot(piecewise.irf_unbounded.h.epsz,'k-','LineWidth',2); hold on;
plot(piecewise.irf.h.epsz,'r-','LineWidth',2); 
title('Hours')
subplot(1,3,3)
plot(piecewise.irf_unbounded.b.epsz,'k-','LineWidth',2); hold on;
plot(piecewise.irf.b.epsz,'r-','LineWidth',2); 
title('Bonds')

% Metrics
moments(1,1) = mean(piecewise_sim.sim.c(burn_in+1:end));
moments(2,1) = mean(piecewise_sim.sim.h(burn_in+1:end));
moments(3,1) = mean(piecewise_sim.sim.b(burn_in+1:end)) /  mean(piecewise_sim.sim.c(burn_in+1:end));
moments(1,2) = std(piecewise_sim.sim.c(burn_in+1:end)) /  mean(piecewise_sim.sim.c(burn_in+1:end));
moments(2,2) = std(piecewise_sim.sim.h(burn_in+1:end)) /  mean(piecewise_sim.sim.h(burn_in+1:end));
moments(3,2) = std(piecewise_sim.sim.b(burn_in+1:end)) /  mean(piecewise_sim.sim.c(burn_in+1:end));
moments(1,3) = skewness(piecewise_sim.sim.c(burn_in+1:end));
moments(2,3) = skewness(piecewise_sim.sim.h(burn_in+1:end));
moments(3,3) = skewness(piecewise_sim.sim.b(burn_in+1:end));

disp('**-- Borrowing constraints model --** ')
disp( table( moments(:,1) , moments(:,2) , moments(:,3) , ...
          'VariableNames',{'Mean','StandardDeviation','Skewness'},...
          'RowNames',{'Consumption';'Hours';'Bonds'}) )

con_binding = zeros(size(piecewise_sim.sim.b(burn_in+1:end)));
con_binding(piecewise_sim.sim.b(burn_in+1:end)+0.01<1e-6) = 1;
disp(['Constraint binds in ', num2str(100*mean(con_binding)),'% of periods'])     
      
% Metrics
moments(1,1) = mean(piecewise_sim.sim_unbounded.c(burn_in+1:end));
moments(2,1) = mean(piecewise_sim.sim_unbounded.h(burn_in+1:end));
moments(3,1) = mean(piecewise_sim.sim_unbounded.b(burn_in+1:end)) / mean(piecewise_sim.sim_unbounded.c(burn_in+1:end));
moments(1,2) = std(piecewise_sim.sim_unbounded.c(burn_in+1:end)) / mean(piecewise_sim.sim_unbounded.c(burn_in+1:end));
moments(2,2) = std(piecewise_sim.sim_unbounded.h(burn_in+1:end)) / mean(piecewise_sim.sim_unbounded.h(burn_in+1:end));
moments(3,2) = std(piecewise_sim.sim_unbounded.b(burn_in+1:end)) / mean(piecewise_sim.sim_unbounded.c(burn_in+1:end));
moments(1,3) = skewness(piecewise_sim.sim_unbounded.c(burn_in+1:end));
moments(2,3) = skewness(piecewise_sim.sim_unbounded.h(burn_in+1:end));
moments(3,3) = skewness(piecewise_sim.sim_unbounded.b(burn_in+1:end));

disp('**-- Linear (unbounded) model --** ')
disp( table( moments(:,1) , moments(:,2) , moments(:,3) , ...
          'VariableNames',{'Mean','StandardDeviation','Skewness'},...
          'RowNames',{'Consumption';'Hours';'Bonds'}) )
    
 
%% Clean up
dynare_cleanup;
