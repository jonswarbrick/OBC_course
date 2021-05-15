%% Solves SOE model using regime-switching and dynare, computing transition probabilities
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
% Burn-in and get random sample
b_samp = sort(randsample(b(burn_in+1:T-1),100)');
z_samp = sort(randsample(z(burn_in+1:T-1),100)');
% Split out regime 1 and regime 2:
b_samp_1 = b_samp(b_samp>b_limit);
z_samp_1 = z_samp(b_samp>b_limit);
z_samp_2 = z_samp(b_samp==b_limit);

p_11 = 1 - normcdf( ( b_limit - b1_cons - g_b1b * ( b_samp_1 - b1_cons ) - g_b1z * ( z_samp_1 - z_cons ) ) / g_b1u );
p_22 = normcdf( ( - mu2_cons - g_mu2z * ( z_samp_2 - z_cons ) ) ./ g_mu2u );
    
h = figure;
subplot(1,2,1)
plot( b_samp_1 , (1-p_11)  , 'k-' , 'LineWidth' , 1.5 );
ylabel('p(1,2)')
xlabel('b(t)')
subplot(1,2,2)
plot( z_samp_1 , (1-p_11)  , 'k-' , 'LineWidth' , 1.5 );
ylabel('p(1,2)')
xlabel('z(t)') 
h.Position = [ 100 , 100 , 800 , 400 ];
sgtitle('Model consistant transition probabilities, p(1,2)')
%print(h,'RS_p11','-r300','-depsc')

h = figure;
plot( z_samp_2 , (1-p_22)  , 'k-' , 'LineWidth' , 1.5 );
ylabel('p(2,1)')
xlabel('z(t)')
h.Position = [ 100 , 100 , 800 , 400 ];
title('Model consistant transition probabilities, p(2,1)')
%print(h,'RS_p22','-r300','-depsc')

% Show RISE model probabilities
psi1 = 500;
psi2 = -200;
theta_1_2 = 2;
theta_2_1 = 2;
p_1_2 = theta_1_2./(theta_1_2+exp(psi1.*(b_samp_1-b_limit)));
p_2_1 = theta_2_1./(theta_2_1+exp(psi2.*z_samp_2));

h = figure;
subplot(1,2,1)
plot( b_samp_1 , (1-p_11)  , 'k--' , 'LineWidth' , 1 ); hold on;
plot( b_samp_1 , p_1_2  , 'r-' , 'LineWidth' , 1.5 );
ylabel('p(1,2)')
xlabel('b(t)')
subplot(1,2,2)
plot( z_samp_2 , (1-p_22)  , 'k--' , 'LineWidth' , 1 ); hold on;
plot( z_samp_2 , p_2_1  , 'r-' , 'LineWidth' , 1.5 );
ylabel('p(2,1)')
xlabel('z(t)') 
legend({'Simulation implied' , 'User defined'})
h.Position = [ 100 , 100 , 800 , 400 ];
sgtitle('User-defined transition probabilities')

%%
dynare_cleanup( );
