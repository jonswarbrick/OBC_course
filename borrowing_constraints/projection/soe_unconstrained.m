%% Solves SOE model using projection
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

disp('Initialisation...')

% Parameters
betta = 0.99;
r = 1/betta;
delta = 0.01;
rho = 0.9;
chi = 0.5;
sigma = 0.01;
y_bar = 1;
b_limit = -0.01;
 
% Projection method settings
Nb = 80;
Nz = 50;
min_b = -1;
max_b = 1;
min_z = -0.1;
max_z = 0.1;
err_tol = 1e-8;
% Nodes and weights for numerical integration
q_pts = 15;
% Grids for the state
b = linspace( min_b , max_b , Nb ); 
z = linspace( min_z , max_z , Nz );
[ b_mesh , z_mesh ] = ndgrid( b , z );
% Simulation settings
time_horizon = 300000;
burn_in = 1000;
plot_start = 5000;
plot_end = 5500;
% IRF simulation settings
irf_periods = 60;
irf_drop = 100;
replic = 100;
IRF_shockscale = 2;

% For plots
mid_z = 25;
low_z = 5;
high_z= 45;

% Initial guess of policy function
b_p = b_mesh;

% Compute expected marginal utility
params = [ rho , sigma , chi , r , q_pts , Nz , Nb ];
muc_p = ex_muc( b_mesh , z_mesh , b_p , params );

% Compute Euler error at each point and sum of square residuals to test for
% convergence
error = muc( exp( z_mesh ) + r * b_mesh - b_p ) - muc_p + delta * b_p;
ssr = sum( sum( error.^2 ) );

disp('Solving policy function...')

% The loop
while ssr>err_tol
    
    % Calculate b_p from RHS of Euler equation
    b_p = exp( z_mesh ) + r * b_mesh -  ( 1 + chi ) ./ ( muc_p - delta * b_p ); 
    muc_p = ex_muc( b_mesh , z_mesh , b_p , params );
    error = muc( ( exp( z_mesh ) + r * b_mesh - b_p ) / (1 + chi) ) - muc_p + delta * b_p;
    ssr = sum( sum( error.^2 ) );
    
end


%% Simulation
disp('Performing time-series simulation...')

% Simulation
rng(147);
eps = randn( 1, time_horizon );
prod = zeros( 1 , time_horizon );
bonds = zeros( 1 , time_horizon );
cons = ones( 1 , time_horizon ) / (1+chi);
for t=2:time_horizon
    prod(t) = rho * prod(t-1) + sigma * eps(t);
    bonds(t) =  interp2( z_mesh , b_mesh , b_p , prod(t) , bonds(t-1) , 'spline' );
end
cons(2:end) = ( exp( prod(2:end) ) + r .* bonds(1:end-1) - bonds(2:end) ) ./ (1 + chi);
hours = 1 - chi .* cons ./ exp( prod );

cons = cons(burn_in+1:end);
hours = hours(burn_in+1:end);
prod = prod(burn_in+1:end);
bonds = bonds(burn_in+1:end);

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

disp( table( moments(:,1) , moments(:,2) , moments(:,3) , ...
          'VariableNames',{'Mean','StandardDeviation','Skewness'},...
          'RowNames',{'Consumption';'Hours';'Bonds'}) )
      
con_binding = zeros(size(bonds));
con_binding(bonds-b_limit<1e-6) = 1;
disp(['Borrowing below constraint  in ', num2str(100*mean(con_binding)),'% of periods'])
      
disp('Simulating average impulse response functions...')
% IRF simulation
irf_horizon = irf_periods+irf_drop;
% Draw shocks
eps_mat = randn( replic, irf_horizon );
% Preallocation
prod_mat = zeros( replic , irf_horizon );
bonds_mat = zeros( replic , irf_horizon );
hours_mat = zeros( replic , irf_horizon );
cons_mat = ones( replic , irf_horizon ) / (1+chi);

% Repeat to get average IRFs
for ii=1:replic
    % Preallocation
    prod_1 = zeros( 1 , irf_horizon );
    bonds_1 = zeros( 1 , irf_horizon );
    cons_1 = ones( 1 , irf_horizon ) / (1+chi);
    prod_2 = prod_1;
    bonds_2 = bonds_1;
    cons_2 = cons_1;
    eps_1 = eps_mat( ii , :);
    eps_2 = eps_1;
    % Add shock of interest
    eps_2( irf_drop+1 ) = eps_2( irf_drop+1 )-IRF_shockscale;
    % Simulation
    for t=2:irf_horizon
        prod_1(t) = rho * prod_1(t-1) + sigma * eps_1(t);
        bonds_1(t) =  interp2( z_mesh , b_mesh , b_p , prod_1(t) , bonds_1(t-1) , 'spline' );
        prod_2(t) = rho * prod_2(t-1) + sigma * eps_2(t);
        bonds_2(t) =  interp2( z_mesh , b_mesh , b_p , prod_2(t) , bonds_2(t-1) , 'spline' );
    end
    cons_1(2:end) = ( exp( prod_1(2:end) ) + r .* bonds_1(1:end-1) - bonds_1(2:end) ) ./ (1 + chi);
    hours_1 = 1 - chi .* cons_1 ./ exp( prod_1 );
    cons_2(2:end) = ( exp( prod_2(2:end) ) + r .* bonds_2(1:end-1) - bonds_2(2:end) ) ./ (1 + chi);
    hours_2 = 1 - chi .* cons_2 ./ exp( prod_2 );
    
    prod_mat( ii , : ) = prod_2 - prod_1;
    bonds_mat( ii , : ) = bonds_2 - bonds_1;
    cons_mat( ii , : ) = cons_2 - cons_1;
    hours_mat( ii , : ) = hours_2 - hours_1;
end

if exist('results/projection_irfs.mat','file'); load('results/projection_irfs.mat'); end

projection.irf_unbounded.c.epsz = mean(cons_mat(:,irf_drop+1:end))/ mean(cons);  % as percent of average consumption
projection.irf_unbounded.b.epsz = mean(bonds_mat(:,irf_drop+1:end))/mean(cons);  % relative to average consumption
projection.irf_unbounded.h.epsz = mean(hours_mat(:,irf_drop+1:end))/ mean(hours); % relative to average hours 

save('results/projection_irfs.mat','projection')

%% Plot results
disp('Plotting...')

figure;
subplot( 2 , 2 , 1)
plot( prod( plot_start:plot_end ) ); title('productivity')
subplot( 2 , 2 , 2)
plot( cons( plot_start:plot_end ) ); title('consumption')
subplot( 2 , 2 , 3)
plot( hours( plot_start:plot_end ) ); title('hours')
subplot( 2 , 2 , 4)
plot( bonds( plot_start:plot_end ) ); title('bonds')
sgtitle('Simulations')

figure;
plot( b , b_p( : , low_z ) ); hold on;
plot( b , b_p( : , mid_z ) )
plot( b , b_p( : , high_z ) )
legend('Low productivity','Middle Productivity','High Productivity')
title('Policy function')
xlabel( 'b' ) 
ylabel( 'bp')

figure;
subplot( 2 , 2 , 1)
plot( projection.irf_unbounded.c.epsz  ); title('consumption')
subplot( 2 , 2 , 2)
plot( projection.irf_unbounded.h.epsz ); title('hours')
subplot( 2 , 2 , 3)
plot( projection.irf_unbounded.b.epsz ); title('bonds')
sgtitle('IRF to technology shock')

disp('All done!')
    
%% Local functions

function marginal_utility = muc( C )

    marginal_utility = 1 ./ C;

end

function muc_p = ex_muc( b_mesh , z_mesh , b_p , params )

    rho =   params( 1 );
    sigma = params( 2 ); 
    chi =   params( 3 ); 
    r =     params( 4 ); 
    q_pts = params( 5 );  

    muc_p  = zeros( size( b_mesh ) );
    [q_n,q_w]= hernodes(q_pts);
    eps_p = sqrt(2)*sigma*q_n;


    for i_q = 1:q_pts

        z_p = rho * z_mesh + eps_p( i_q );
        b_p_p = interp2( z_mesh , b_mesh , b_p , z_p , b_p , 'spline' );
        muc_p = muc_p +  q_w(i_q) * muc( ( exp( z_p ) + r * b_p - b_p_p ) / ( 1 + chi ) );

    end
    muc_p = muc_p ./ sqrt(pi);

end