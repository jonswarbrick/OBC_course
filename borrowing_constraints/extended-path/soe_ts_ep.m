%% Computes extended-path simulation for SOE model with borrowing constraint
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
ss.u = log( ss.c ) + p.chi * log( 1 - ss.h ) - p.delta * ss.b^2;
ss.v = ss.u / ( 1 - p.betta );


% Solution and simulation settings
T = 1000; % overall including burn in and out
burn_in = 100;
burn_out = 100;
crit = 1e-12;
err = 1;

% Draw shocks
z = zeros(T,1);
rng(147)
epsz = randn( T, 1);
epsz(T-burn_out+1:end) = zeros( burn_out , 1 ); % Set burn in period to zero
z(1) = p.sigma * epsz(1);
for t=2:T
    z(t) = p.rho * z(t-1) + p.sigma * epsz(t);
end

% Extended path simulation
disp('Computing extended-path simulation')
options = optimoptions('fsolve','Display','none');
b0 = ss.b * ones( T , 1 ); % Initial values
b = b0; % Initialization

% Period-by-period computation
for t=2:T
    disp(['Computation for expectated path at time ',num2str(t),' of ',num2str(T)])
    
    b0 = b(t:end);
    
    fun = @(bi) model_f( bi , b , z , ss , p , t ); % The model (obc)
    [bi , ~ , flag , ~] = fsolve( fun , b0 , options );
    
    b = [ b(1:t-1) ; bi ];
    
    if flag<1
        disp('**--- ERROR ---**')
    end
end
% Consumption and hours
ep.b = b;
ep.c = ( exp( z ) + p.r .* [ ss.b ; ep.b(1:end-1) ] -  ep.b ) ./ ( 1 + p.chi );
ep.h = 1 - p.chi .* ep.c ./ exp(z);

ep.b = ep.b(burn_in+1:end-burn_out);
ep.c = ep.c(burn_in+1:end-burn_out);
ep.h = ep.h(burn_in+1:end-burn_out);


% Perfect foresight simulation
disp('Computing perfect foresight simulation')
fun = @(b) model_f_pf( b , z , ss , p ); % The model (obc)
disp('Computing perfect-foresight simulation')
options = optimoptions('fsolve','Display','iter');
b0 = ss.b * ones( T , 1 ); % Initial values 
[pf.b , ~ , flag , ~] = fsolve( fun , b0 , options );
if flag>0
    disp('Solved')
    
    % Consumption and hours
    pf.c = ( exp( z ) + p.r .* [ ss.b ; pf.b(1:end-1) ] -  pf.b ) ./ ( 1 + p.chi );
    pf.h = 1 - p.chi .* pf.c ./ exp(z);
    
    pf.b = pf.b(burn_in+1:end-burn_out);
    pf.c = pf.c(burn_in+1:end-burn_out);
    pf.h = pf.h(burn_in+1:end-burn_out);

else
    disp('Not solved')
end

z = z(burn_in+1:end-burn_out);

%% Plotting
disp('Plotting...')

h = figure;
subplot( 2 , 2 , 1)
plot( z , 'k-' , 'LineWidth' , 1 ); 
title('z')
subplot( 2 , 2 , 2)
plot( ep.c , 'k-' , 'LineWidth' , 1 ); hold on;
plot( pf.c , 'Color' , [.6 0 .2] , 'LineWidth' , 1 , 'LineStyle' , '-' );
title('c')
legend({'Extended path','Perfect foresight'})
subplot( 2 , 2 , 3)
plot( ep.h , 'k-' , 'LineWidth' , 1 ); hold on;
plot( pf.h , 'Color' , [.6 0 .2] , 'LineWidth' , 1 , 'LineStyle' , '-' );
title('h')
subplot( 2 , 2 , 4)
plot( ep.b , 'k-' , 'LineWidth' , 1 ); hold on;
plot( pf.b , 'Color' , [.6 0 .2] , 'LineWidth' , 1 , 'LineStyle' , '-' );
title('b')
sgtitle('Time-series simulations')
h.Position = [ 100 , 100 , 800 , 400 ];


plot_t = 600;
% Plot expected shocks
for t=1:10

    s=plot_t+t-1;
    
    zt = z;
    % Compute path for z if no more shocks
    for j=s+1:length(zt)
        zt(j) = p.rho * zt(j-1);
    end
    expz.(['z',num2str(t)]) = zt;
    
end

h = figure;
plot(0:50,expz.z1(plot_t:plot_t+50) ); hold on;
plot(0:50,expz.z2(plot_t:plot_t+50) )
plot(0:50,expz.z3(plot_t:plot_t+50) )
plot(0:50,expz.z4(plot_t:plot_t+50) )
plot(0:50,expz.z5(plot_t:plot_t+50) )
plot(0:50,expz.z6(plot_t:plot_t+50) )
plot(0:50,expz.z7(plot_t:plot_t+50) )
plot(0:50,expz.z8(plot_t:plot_t+50) )
plot(0:50,expz.z9(plot_t:plot_t+50) )
plot(0:50,expz.z10(plot_t:plot_t+50) )
plot(0:50,z(plot_t:plot_t+50) , 'k-' , 'LineWidth' , 1 )
h.Position = [ 100 , 100 , 800 , 400 ];
sgtitle('Shock Expectations')
ylabel('z')

%% Local functions
function F = model_f( bi , b , z , ss , p , t ) 

    remain_t = length(bi);
    lag_b = [ b(t-1) ; bi(1:end-1) ];
    lead_b = [ bi(2:end) ; ss.b ];
    zi = zeros( size( bi ) );
    zi(1) = z(t);
    % Compute path for z if no more shocks
    for s=2:remain_t
        zi(s) = p.rho * zi(s-1);
    end
    lead_z = [ zi(2:end) ; 0 ];
    c = ( exp( zi ) + p.r .* lag_b - bi ) ./ (1 + p.chi);
    lead_c = ( exp( lead_z ) + p.r .* bi - lead_b ) ./ (1 + p.chi);

    mu = (1./c) - p.betta * p.r * (1./lead_c) + p.delta .* bi;
   
    F = min( mu , bi-p.b_limit );

end

function F = model_f_pf( b , z , ss , p ) 

    lag_b = [ ss.b ; b(1:end-1) ];
    c = ( exp( z ) + p.r .* lag_b - b ) ./ (1 + p.chi);
    lead_c = [ c(2:end) ; ss.c ];

    mu = (1./c) - p.betta * p.r * (1./lead_c) + p.delta .* b;
   
    F = min( mu , b-p.b_limit );

end
