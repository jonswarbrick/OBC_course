%% Solves SOE model with a borrowing constraint using a piecewise linear approach
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

%% Run dynare to get 1st order policy function
dynare soe_unconstrained

b_ss = oo_.dr.ys(strmatch('b',M_.endo_names,'exact'));
z_ss = oo_.dr.ys(strmatch('z',M_.endo_names,'exact'));
gbb = oo_.dr.ghx((oo_.dr.order_var==strmatch('b',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('b',M_.endo_names,'exact')) );
gbz = oo_.dr.ghx((oo_.dr.order_var==strmatch('b',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );
gzb = oo_.dr.ghx((oo_.dr.order_var==strmatch('z',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('b',M_.endo_names,'exact')) );
gzz = oo_.dr.ghx((oo_.dr.order_var==strmatch('z',M_.endo_names,'exact')) , (oo_.dr.state_var==strmatch('z',M_.endo_names,'exact')) );

gbu = oo_.dr.ghu((oo_.dr.order_var==strmatch('b',M_.endo_names,'exact')));
gzu = oo_.dr.ghu((oo_.dr.order_var==strmatch('z',M_.endo_names,'exact')));


% Parameters
p.betta = 0.99;
p.r = 1/p.betta;
p.delta = 0.01;
p.rho = 0.9;
p.chi = 0.5;
p.sigma = 0.01;
p.y_bar = 1;
p.b_limit = -0.01;

ss.z = 0;
ss.b = 0;
ss.c = 1 / ( 1 + p.chi );
ss.h = ss.c; 

% Solution and simulation settings
H = 200; 

disp('Begin piecewise-linear simulation')
% Draw shock
b = b_ss * ones(H,1);
blin = b_ss * ones(H,1);
z = zeros(H,1);
epsz = zeros(H,1);
epsz(2) = - 2; %shock in period 2
% Presimulate (exogenous) z and linear b
for t=2:H
    z(t) = z_ss + gzz * z(t-1) + gzu * epsz(t);
    blin(t) = b_ss + gbb * blin(t-1) + gbz * z(t-1) + gbu * epsz(t);
end
t = 2;
chk_1 = 1;
T_guess = 5;
while chk_1==1
    
    b(t) = b_ss + gbb * b(t-1) + gbz * z(t-1) + gbu * epsz(t);
    if b(t)<p.b_limit
        % Guess T to return to steady state
        disp(['Constraint violating, trying first guess T = ',num2str(T_guess)]);
        chk_2 = 1;
        T = T_guess;
        while chk_2==1
            b_T = b_ss + gbb * p.b_limit + gbz * z(t-1+T) + gbu * epsz(t+T);
            lag_b_T = b_ss + gbb * p.b_limit + gbz * z(t-1+T-1) + gbu * epsz(t+T-1);
            if b_T<p.b_limit
                % If OBC still binding, increase horizon
                disp(['Increasing T, trying T = ',num2str(T+1)]);
                T = T+1;
            elseif lag_b_T>p.b_limit
                % Check if the horizon should be shorter
                disp(['Decreasing T, trying T = ',num2str(T-1)]);
                T = T-1;
            else
                b(t+T) = b_ss + gbb * p.b_limit + gbz * z(t-1+T) + gbu * epsz(t+T);
                for s=1:T
                    b(t+T-s) = p.b_limit;
                end
                t=t+T;
                chk_2 = 0;
            end
        end
    end
    if t==H
        chk_1 = 0;
    else
        t=t+1;
    end
end

c = ( exp( z ) + p.r .* [ ss.b ; b(1:end-1) ] -  b ) ./ ( 1 + p.chi );
h = 1 - p.chi .* c ./ exp(z);
clin = ( exp( z ) + p.r .* [ ss.b ; blin(1:end-1) ] -  blin ) ./ ( 1 + p.chi );
hlin = 1 - p.chi .* clin ./ exp(z);

figure;
subplot( 2 , 2 , 1)
plot( clin(2:60) , 'k-' , 'LineWidth' , 1.5 ); hold on;
plot( c(2:60) , 'r-' , 'LineWidth' , 1.5 ); 
title('c')
subplot( 2 , 2 , 2)
plot( hlin(2:60) , 'k-' , 'LineWidth' , 1.5 );  hold on;
plot( h(2:60) , 'r-' , 'LineWidth' , 1.5 ); 
title('h')
legend({'Linear','OBC'})
subplot( 2 , 2 , 3)
plot( blin(2:60) , 'k-' , 'LineWidth' , 1.5 );  hold on;
plot( b(2:60) , 'r-' , 'LineWidth' , 1.5 ); 
title('b')

%% Clean up
dynare_cleanup;