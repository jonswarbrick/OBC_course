%% Computes perfect-foresight simulation for SOE model with borrowing constraint
% For the course "Occasionally Binding Constraints in Macro"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

disp('Initialisation...')
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
ss.u = log( ss.c ) + p.chi * log( 1 - ss.h ) - 2 * p.delta * ss.b^2;
ss.v = ss.u / ( 1 - p.betta );


% Solution and simulation settings
T = 1000; % overall including burn in and out
burn_in = 100;
burn_out = 100;
crit = 1e-12;
err = 1;

% Draw shocks
z = zeros(T,1);
rng(147); 
epsz = randn( T, 1);
epsz(1:burn_in) = zeros( burn_in , 1 ); % Set burn in period to zero
epsz(T-burn_out+1:end) = zeros( burn_out , 1 ); % Set burn in period to zero
z(1) = p.sigma * epsz(1);
for t=2:T
    z(t) = p.rho * z(t-1) + p.sigma * epsz(t);
end

fun = @(b) model_f( b , z , ss , p ); % The model (obc)
fun_unc = @(b) model_f_unc( b , z , ss , p ); % The model (unconstrained) 

disp('Computing path of endogenous state for OBC model')
options = optimoptions('fsolve','Display','iter');
b0 = ss.b * ones( T , 1 ); % Initial values 
[obc.b , ~ , flag , ~] = fsolve( fun , b0 , options );
if flag>0
    disp('Solved')
    
    % Compute consumption and hours
    obc.c = ( exp( z ) + p.r .* [ ss.b ; obc.b(1:end-1) ] -  obc.b ) ./ ( 1 + p.chi );
    obc.h = 1 - p.chi .* obc.c ./ exp(z);
    
    obc.b = obc.b(burn_in+1:end-burn_out);
    obc.c = obc.c(burn_in+1:end-burn_out);
    obc.h = obc.h(burn_in+1:end-burn_out);

else
    disp('Problem finding solution to OBC model')
end


disp('Computing path of endogenous state for unconstrained model')
options = optimoptions('fsolve','Display','iter');
b0 = ss.b * ones( T , 1 ); % Initial values 
[unc.b , ~ , flag , ~] = fsolve( fun_unc , b0 , options );
if flag>0
    disp('Solved')
    
    % Compute consumption and hours
    unc.c = ( exp( z ) + p.r .* [ ss.b ; unc.b(1:end-1) ] -  unc.b ) ./ ( 1 + p.chi );
    unc.h = 1 - p.chi .* unc.c ./ exp(z);
    
    unc.b = unc.b(burn_in+1:end-burn_out);
    unc.c = unc.c(burn_in+1:end-burn_out);
    unc.h = unc.h(burn_in+1:end-burn_out);
else
    disp('Problem finding solution to OBC model')
end

z = z(burn_in+1:end-burn_out);


%% Plotting
disp('Plotting...')

h = figure;
subplot( 2 , 2 , 1)
plot( z , 'k-' , 'LineWidth' , 1 ); hold on;
plot( z , 'Color' , [.6 0 .2] , 'LineWidth' , 1 );
title('z')
legend({'OBC','No OBC'})
subplot( 2 , 2 , 2)
plot( obc.c , 'k-' , 'LineWidth' , 1 ); hold on;
plot( unc.c , 'Color' , [.6 0 .2] , 'LineWidth' , 1 );
title('c')
subplot( 2 , 2 , 3)
plot( obc.h , 'k-' , 'LineWidth' , 1 ); hold on;
plot( unc.h , 'Color' , [.6 0 .2] , 'LineWidth' , 1 );
title('h')
subplot( 2 , 2 , 4)
plot( obc.b , 'k-' , 'LineWidth' , 1 ); hold on;
plot( unc.b , 'Color' , [.6 0 .2] , 'LineWidth' , 1 );
title('b')
sgtitle('Time-series simulations')
h.Position = [ 100 , 100 , 800 , 400 ];


%% Local functions
function F = model_f( b , z , ss , p ) 

    lag_b = [ ss.b ; b(1:end-1) ];
    c = ( exp( z ) + p.r .* lag_b - b ) ./ (1 + p.chi);
    lead_c = [ c(2:end) ; ss.c ];

    mu = (1./c) - p.betta * p.r * (1./lead_c) + 2.* p.delta .* b;
   
    F = min( mu , b-p.b_limit );

end

function F = model_f_unc( b , z , ss , p ) 

    lag_b = [ ss.b ; b(1:end-1) ];
    c = ( exp( z ) + p.r .* lag_b - b ) ./ (1 + p.chi);
    lead_c = [ c(2:end) ; ss.c ];

    mu = (1./c) - p.betta * p.r * (1./lead_c) + 2.* p.delta .* b;
   
    F = mu;

end

