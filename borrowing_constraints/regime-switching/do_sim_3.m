%% Solves SOE model using regime-switching and dynare 1-step iteration on transition probabilities
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

% Run dynare to get model structs and initial policy functions
dynare soe_regime_switching

%% Solve transition probabilities

disp('//**---------------------------**//')
disp('Retreiving policy functions')
% Get policy function coefficients
b1_cons = oo_.dr.ys(strmatch('b1',M_.endo_names,'exact'));
g_b1b = oo_.dr.ghx((oo_.dr.order_var==strmatch('b1',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('b1',M_.endo_names,'exact')) );
g_b1z = oo_.dr.ghx((oo_.dr.order_var==strmatch('b1',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );
g_b1u = oo_.dr.ghu((oo_.dr.order_var==strmatch('b1',M_.endo_names,'exact')));
z_cons = oo_.dr.ys(strmatch('z',M_.endo_names,'exact'));
g_zz = oo_.dr.ghx((oo_.dr.order_var==strmatch('z',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );
g_zu = oo_.dr.ghu((oo_.dr.order_var==strmatch('z',M_.endo_names,'exact')));
mu2_cons = oo_.dr.ys(strmatch('mu2',M_.endo_names,'exact'));
g_mu2z = oo_.dr.ghx((oo_.dr.order_var==strmatch('mu2',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );
g_mu2u = oo_.dr.ghu((oo_.dr.order_var==strmatch('mu2',M_.endo_names,'exact')));


% Simulation settings
% Simulate model over 10000 periods, with 100 period burn-in
T = 300000;
burn_in = 1000;
b_limit = M_.params(strmatch('b_limit',M_.param_names,'exact'));
R = M_.params(strmatch('R',M_.param_names,'exact'));
chi = M_.params(strmatch('chi',M_.param_names,'exact'));

% Perform simulation to recover transition probabilities
b = b1_cons * ones( 1 , T );
b1 = b;
mu2 = mu2_cons * ones( 1 , T );
z = zeros( 1 , T );
epsz = randn( 1 , T );
flag = 1; % Indicate which regime
for t=2:T
    z(t) = z_cons + g_zz * z(t-1) + g_zu * epsz(t);    
    b1(t) = b1_cons + g_b1b * ( b(t-1) - b1_cons ) + g_b1z * z(t-1) + g_b1u * epsz(t);
    mu2(t) = mu2_cons + g_mu2z * z(t-1) + g_mu2u * epsz(t);

    switch(flag)
        case 1
            if b1(t)<b_limit
                b1(t) = b_limit;
                b(t) = b_limit;
                flag = 2;
            else
                b(t) = b1(t);
            end
        case 2
            if mu2(t)<0
                b(t) = b1(t);
                flag = 1;
            else
                b1(t) = b_limit;
                b(t) = b_limit;
            end
    end
end
mean_z2 = mean(z(b==b_limit));

% Regime 1 transition probabilities, linear coefficients
% Compute around deterministic steady state
theta11_0_ = 1 - normcdf( ( b_limit - b1_cons ) / g_b1u );
theta11_1_ = diff( 1 - normcdf( ( b_limit - b1_cons - g_b1b * ( [ -0.001 0.001 ] ) ) ./ g_b1u ) ) / diff( [ -0.001 0.001 ] );
theta11_2_ = diff( 1 - normcdf( ( b_limit - b1_cons - g_b1z * ( [ -0.001 0.001 ] ) ) ./ g_b1u ) ) / diff( [ -0.001 0.001 ] );

% Regime 2 transition probabilities, linear coefficients
% Compute around conditional (b=b_limit) ergodic mean
theta22_0_ = normcdf( ( - mu2_cons - g_mu2z * ( mean( mean_z2 ) - z_cons ) ) ./ g_mu2u );
theta22_2_ = diff( normcdf(  - ( mu2_cons + g_mu2z * ( mean_z2 + [ -0.001 0.001 ] ) ) ./ g_mu2u ) ) / diff( [ -0.001 0.001 ] );


fprintf('\nRe-solving model and retrieving updated policy function \n')
set_param_value( 'theta11_0' , theta11_0_ );
set_param_value( 'theta11_1' , theta11_1_ );
set_param_value( 'theta11_2' , theta11_2_ );
set_param_value( 'theta22_0' , theta22_0_ );
set_param_value( 'theta22_2' , theta22_2_ );
[oo_.dr ,~,M_,~,oo_] = resol(0,M_,options_ ,oo_);
    
% Get updated policy function coefficients
b1_cons = oo_.dr.ys(strmatch('b1',M_.endo_names,'exact'));
g_b1b = oo_.dr.ghx((oo_.dr.order_var==strmatch('b1',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('b1',M_.endo_names,'exact')) );
g_b1z = oo_.dr.ghx((oo_.dr.order_var==strmatch('b1',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );
g_b1u = oo_.dr.ghu((oo_.dr.order_var==strmatch('b1',M_.endo_names,'exact')));
z_cons = oo_.dr.ys(strmatch('z',M_.endo_names,'exact'));
g_zz = oo_.dr.ghx((oo_.dr.order_var==strmatch('z',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );
g_zu = oo_.dr.ghu((oo_.dr.order_var==strmatch('z',M_.endo_names,'exact')));
mu2_cons = oo_.dr.ys(strmatch('mu2',M_.endo_names,'exact'));
g_mu2z = oo_.dr.ghx((oo_.dr.order_var==strmatch('mu2',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );
g_mu2u = oo_.dr.ghu((oo_.dr.order_var==strmatch('mu2',M_.endo_names,'exact')));
   
% New simulation 
fprintf('\nPerform simulation with updated linear transition probabilities \n')
b = b1_cons * ones( 1 , T );
b1 = b;
mu2 = mu2_cons * ones( 1 , T );
z = zeros( 1 , T );
epsz = randn( 1 , T );
flag = 1; % Indicate which regime
for t=2:T
    z(t) = z_cons + g_zz * z(t-1) + g_zu * epsz(t);    
    b1(t) = b1_cons + g_b1b * ( b(t-1) - b1_cons ) + g_b1z * z(t-1) + g_b1u * epsz(t);
    mu2(t) = mu2_cons + g_mu2z * z(t-1) + g_mu2u * epsz(t);

    switch(flag)
        case 1
            if b1(t)<b_limit
                b1(t) = b_limit;
                b(t) = b_limit;
                flag = 2;
            else
                b(t) = b1(t);
            end
        case 2
            if mu2(t)<0
                b(t) = b1(t);
                flag = 1;
            else
                b1(t) = b_limit;
                b(t) = b_limit;
            end
    end
end


disp('Complete!')

c = ( exp( z ) +  R .* [ b1_cons , b(1:end-1) ] -  b ) ./ ( 1 + chi );
h = 1 - chi .* c ./ exp(z);

cons = c(burn_in+1:end);
hours = h(burn_in+1:end);
prod = z(burn_in+1:end);
bonds = b(burn_in+1:end);

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

%% Clean up
dynare_cleanup;