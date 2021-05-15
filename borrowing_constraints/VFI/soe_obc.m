%% Solves SOE model with a borrowing constraint using Value Function Iteration
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')
disp('Initialisation...')

% Model parameters
p.betta = 0.99;
p.r = 1/p.betta;
p.delta = 0.01;
p.rho = 0.9;
p.chi = 0.5;
p.sigma = 0.01;
p.b_limit = -0.01;

% Solution settings
b_pts = 40;
b_p_pts = 1000;
z_pts = 40;
min_b = p.b_limit;
max_b = 1;
min_z = -0.1; 
max_z = 0.1;
% Nodes and weights for numerical integration
q_pts = 15;
err_tol = 1e-7;

% Simulation settings
time_horizon = 10000; % 300000; % for time-series simulation
burn_in = 1000;
irf_periods = 60;
irf_drop = 100;
replic = 100;
IRF_shockscale = 2;

% Plotting options
plot_start = 5000;
plot_end = 5500;
low_z = 5;
mid_z = 20;
high_z= 35;
color_1 = [0, 0.4470, 0.7410];
color_2 = [0, 0.5, 0];
color_3 = [153 0 51] / 255;

% Initializations
%b = linspace(min_b , max_b , b_pts);
%b_p = linspace(min_b , max_b , b_p_pts);
b = exp( linspace( log(min_b+0.02) , log(max_b+0.02) , b_pts ) ) - 0.02; % concentrate points for small b
b_p = exp( linspace( log(min_b+0.02) , log(max_b+0.02) , b_p_pts ) ) - 0.02;  % concentrate points for small b
z = linspace( min_z , max_z , z_pts ); % Could alternatively use tauchen method to choose z points 
[b_mesh , z_mesh , b_p_mesh ] = ndgrid( b , z , b_p ); % Creates grids for maximization step
[b_grid , z_grid ] = ndgrid( b , z ); % Create grids for current value

steady.z = 0;
steady.b = 0;
steady.c = 1 / ( 1 + p.chi );
steady.h = steady.c; 
steady.u = utility( steady.b , steady.b , steady.z , p );
steady.v = steady.u / ( 1 - p.betta );

V = steady.v * ones( b_pts , z_pts ); % Value function first guess
crit = 1000; % Initial value for error criterion
iter = 0; % Initial  iteration number
start_time = cputime; % For timing

disp('Solving value function...')

while crit>err_tol
    
    % cont_V is expected continuation value given current productivity z and
    % choice of b. This integrates over next period z.
    [q_n,q_w]= hernodes(q_pts);
    p_i = ( q_w ./ sqrt(pi) );
    eps_p = sqrt(2)*p.sigma*q_n;
    cont_V = zeros( b_pts , z_pts , b_p_pts );
    for i = 1:q_pts
        z_p_grid = p.rho .* z_mesh + eps_p( i );
        cont_V = cont_V + p.betta * p_i( i ) * V_fun( z_grid , b_grid , V , z_p_grid , b_p_mesh );
    end
   
    % Get maximum continuation value conditional on current value
    % Grid for z and current and future b
    % tv1 is the value of V given optimal choice of b
    % tdecis is the grid indices of optimal choice of b 
    [V_new,b_opt_ind] = max( utility( b_mesh , b_p_mesh , z_mesh , p ) + cont_V , [] , 3 );
    
    % Check for convergence (max relative error)
    crit = max( max( abs( ( V_new - V ) ./ V_new ) ) );
    
    % Update V grid
    V = V_new;
    
    % Report iteration details
    iter = iter+1;
    fprintf('iter = %g ; error criterion = %e\n' , iter , crit);
    
end

disp('  ');
fprintf('Computation time = %f seconds \n', cputime-start_time);

% Transform decision indices into bond choice
b_policy = b_p( b_opt_ind );

h = figure;
plot( b' , b' , 'Color' , [0.7 0.7 0.7 ],'HandleVisibility','off' );hold on
plot( b' , b_policy( : , low_z ) , 'Color' , color_1 , 'LineWidth' , 1.5 );  
plot( b' , b_policy( : , mid_z ) , 'Color' , color_2 , 'LineWidth' , 1.5 )
plot( b' , b_policy( : , high_z ) , 'Color' , color_3 , 'LineWidth' , 1.5 )
legend({'Low productivity','Middle Productivity','High Productivity'},'Location','northwest')
title('Policy function')
xlim([-0.01,0.1])
xlabel( 'b' )
ylabel( 'bp')
h.Position = [ 100 , 100 , 800 , 400 ];
print(h,'VFI_policy','-r300','-depsc')


h = figure;
plot( b' , V(:,low_z) , 'Color' , color_1 , 'LineWidth' , 1.5 );  hold on;
plot( b' , V(:,mid_z) , 'Color' , color_2 , 'LineWidth' , 1.5 );  
plot( b' , V(:,high_z)  , 'Color' , color_3 , 'LineWidth' , 1.5 );  
legend({'z_0','z_1','z_3'},'Location','northwest')
title('Value function')
xlabel( 'b' )
ylabel( 'V(b,z)')
xlim([-0.01,0.2])
h.Position = [ 100 , 100 , 800 , 400 ];
print(h,'VFI_value','-r300','-depsc')


%% Simulation
disp('Performing time-series simulation...')

% Simulation
rng(147);
sim.eps = randn( 1, time_horizon );
sim.z = zeros( 1 , time_horizon );
sim.b = zeros( 1 , time_horizon );
sim.c = ones( 1 , time_horizon ) / (1+p.chi);
for t=2:time_horizon
    sim.z(t) = p.rho * sim.z(t-1) + p.sigma * sim.eps(t);
    sim.b(t) =  max( p.b_limit , interp2( repmat( z , b_pts , 1 ) , repmat( b' , 1 , z_pts ) , b_policy , sim.z(t) , sim.b(t-1) , 'spline' ) );
end
sim.c(2:end) = ( exp( sim.z(2:end) ) + p.r .* sim.b(1:end-1) - sim.b(2:end) ) ./ (1 + p.chi);
sim.h = 1 - p.chi .* sim.c ./ exp( sim.z () );

sim.c = sim.c(burn_in+1:end);
sim.h = sim.h(burn_in+1:end);
sim.z = sim.z(burn_in+1:end);
sim.b = sim.b(burn_in+1:end);

% Metrics
moments(1,1) = mean(sim.c);
moments(2,1) = mean(sim.h);
moments(3,1) = mean(sim.b/mean(sim.c));
moments(1,2) = std(sim.c/mean(sim.c));
moments(2,2) = std(sim.h/mean(sim.h));
moments(3,2) = std(sim.b/mean(sim.c));
moments(1,3) = skewness(sim.c);
moments(2,3) = skewness(sim.h);
moments(3,3) = skewness(sim.b);

disp( table( moments(:,1) , moments(:,2) , moments(:,3) , ...
          'VariableNames',{'Mean','StandardDeviation','Skewness'},...
          'RowNames',{'Consumption';'Hours';'Bonds'}) )

  
con_binding = zeros(size(sim.b));
con_binding(sim.b-p.b_limit<1e-8) = 1;
disp(['Constraint binds in ', num2str(100*mean(con_binding)),'% of periods'])

figure;
subplot( 2 , 2 , 1)
plot( sim.z( plot_start:plot_end ) ); title('productivity')
subplot( 2 , 2 , 2)
plot( sim.c( plot_start:plot_end ) ); title('consumption')
subplot( 2 , 2 , 3)
plot( sim.h( plot_start:plot_end ) ); title('hours')
subplot( 2 , 2 , 4)
plot( sim.b( plot_start:plot_end ) ); title('bonds')
sgtitle('Simulations')


%% Functions

function u = utility( b , b_p , z , p )

c = max( 1e-20 , ( exp( z ) + p.r .* b - b_p ) ./ ( 1 + p.chi ) );
h = 1 - p.chi .* c ./ exp( z );
u = log( c ) + p.chi .* log( 1 - h ) - p.delta .* b_p.^2;

end

function V_p = V_fun( z , b , V , z_p , b_p )

V_p = interp2( z , b , V , z_p , b_p , 'spline' );

end