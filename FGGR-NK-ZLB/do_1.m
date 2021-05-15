%% Computes IRF for FVetal2015 NK model using dynareOBC perfect foresight solver,
% printing news shocks to command window
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../functions')
format long

dynareOBC NK_IRF DisplayBoundsSolutionProgress ShockScale=10

%% Plot IRF using solution to the bounds problem
r_steady = log( 1.005 / 0.994 );
q = dynareOBC_.IRFsWithoutBounds.r_epsilon_b' + r_steady;
M = dynareOBC_.MMatrix( 1:40, 1:5 );
y =  [ ...
   0.055238204233548
   0.031473406542942
   0.016560450270148
   0.007452116533236
   0.002075504046452
    ];
figure; plot( 1:40, q, ':k', 1:40, q + M * y, '-k' );

%% Plot the M Matrix
figure;
plot( M( : , 1 ) );
pause;
plot( M( : , 2 ) );
pause;
plot( M( : , 3 ) );
pause;
plot( M( : , 4 ) );
pause;
plot( M( : , 5 ) );
pause;
plot( dynareOBC_.MMatrix );

%% Plot the news shocks and final solution
figure;
plot( 1:40 , q , 1:40 , M * diag( y ) , 1:40 , q + M * y );

%%
format short
dynare_cleanup( );
