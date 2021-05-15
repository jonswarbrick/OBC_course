%% Solves SOE model with a borrowing constraint using add factors
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../functions')



%% Run dynare to get 1st order policy function with news shocks
dynare soe_addfactors
% With model Y(t) = A*X(t)
% where X(t) includes (t) shocks, news shocks and (t-1) states
A = disp_dr_(oo_.dr,options_.order,[]);

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
TimeToEscapeBounds = 30;


% Preallocate variables N X T
b = ss.b * ones( 1 , H );
z = zeros( 1 , H );
add_factors = zeros( 1 , H );
news_shock =  zeros( 1+TimeToEscapeBounds , H+TimeToEscapeBounds );

% Index in decision rule
z_i = 3;
b_i = 2;
add_factors_i = 4;

% Get state variable index and names
g_ = getStateIndex( TimeToEscapeBounds );
X_list = getStateNames( TimeToEscapeBounds );


disp('Linear simulation')
% Draw shock
epsz = zeros(H,1);
epsz(2) = - 2; %shock in period 2
% Presimulate (exogenous) z and linear b
for t=2:H
    z( t ) = ss.z;
    for ii=1:length(X_list)
        z( t ) = z( t ) + A( ii+1 , z_i ) * eval(char(X_list( ii )));
    end    
    b( t ) = ss.b;
    for ii=1:length(X_list)
        b( t ) = b( t ) + A( ii+1 , b_i ) * eval(char(X_list( ii )));
    end    
end

blin = b;
clin = ( exp( z ) + p.r .* [ ss.b , blin(1:end-1) ] -  blin ) ./ ( 1 + p.chi );
hlin = 1 - p.chi .* clin ./ exp(z);

disp('News shock calculation')
t=1;
while t<H
    t=t+1;
    b( t ) = ss.b;
    for ii=1:length(X_list)
        b( t ) = b( t ) + A( ii+1 , b_i ) * eval(char(X_list( ii )));
    end   
    if b( t )<p.b_limit
        
        % Find news sequence of news shocks
        f = @(y) sim_b( y , b , z , epsz , news_shock , X_list , A , b_i , add_factors_i , t , TimeToEscapeBounds , p );
        nonlincon = @(ys) fcons( ys , b , z , epsz , news_shock , X_list , A , b_i , add_factors_i , t , TimeToEscapeBounds , p );
        y0 = zeros( TimeToEscapeBounds+1 , 1 ); y0(1:18) = [ 0.001839006447162;0.0017683;0.0015714;0.001394323;0.0012348907;0.001091401630000;0.000962261467000;0.000846035320300;0.000741431788270;0.000647288609443;0.000562559748499;0.000486303773649;0.000417673396284;0.000355906056656;0.000300315450990;0.000250283905891;0.000205255515302;0.00004065503976];   
        y = fmincon(f,y0,[],[],[],[],[],[],nonlincon,optimoptions('fmincon','Display','iter','Algorithm','sqp'));
        news_shock( : , t+TimeToEscapeBounds ) = y;  
        b( t ) = p.b_limit;
    end
end
s=3;
add_factor = zeros( size( b ) );
for t=s:s+TimeToEscapeBounds
    for ii=1:length(X_list)
        b( t ) = b( t ) + A( ii+1 , b_i ) * eval(char(X_list( ii )));
        add_factor( t ) = add_factor( t ) + A( ii+1 , add_factors_i ) * eval(char(X_list( ii )));
        
    end
end

c = ( exp( z ) + p.r .* [ ss.b , b(1:end-1) ] -  b ) ./ ( 1 + p.chi );
h= 1 - p.chi .* c ./ exp(z);

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

%% Local functions
function g = getStateIndex( T )

ind=1;
for ii=1:T
    for jj=1:ii
        ind=ind+1;
        g.(['news_',num2str(ii),'_lag',num2str(jj)]) = ind;
    end
end
ind=ind+1;
g.b_lag_1 = ind;
ind=ind+1;
g.z_lag_1 = ind;
ind=ind+1;
g.epsz = ind;
ind=ind+1;
g.news = ind;
for ii=1:T
    ind=ind+1;
    g.(['news_',num2str(ii)]) = ind;
end


end

function X = getStateNames( T )

X = {};
for ii=2:T+1
    for jj=1:ii-1
        X = [ X ; {['news_shock(',num2str(ii),',t+TimeToEscapeBounds-',num2str(jj),')']} ];
    end
end
X = [ X ; [ {'b(t-1)'} ; {'z(t-1)'} ; {'epsz(t)'} ; {'news_shock(1 , t+TimeToEscapeBounds)'}] ];
for ii=2:T+1
    X = [ X ; {['news_shock(',num2str(ii),',t+TimeToEscapeBounds)']} ];
end


end

function F = sim_b( y , b , z , epsz , news_shock , X_list , A , b_i , add_factors_i , s , TimeToEscapeBounds , p )

news_shock( : , s+TimeToEscapeBounds ) = y;

b( s ) = 0;
add_factor = zeros( size( b ) );
for t=s:s+TimeToEscapeBounds
    for ii=1:length(X_list)
        b( t ) = b( t ) + A( ii+1 , b_i ) * eval(char(X_list( ii )));
        add_factor( t ) = add_factor( t ) + A( ii+1 , add_factors_i ) * eval(char(X_list( ii )));
    end
end

F = sum( (( b(s:s+TimeToEscapeBounds)-p.b_limit )'.* y).^2 );

end

function [c,ceq] = fcons( y , b , z , epsz , news_shock , X_list , A , b_i , add_factors_i , s , TimeToEscapeBounds , p )

news_shock( : , s+TimeToEscapeBounds ) = y;

b( s ) = 0;
add_factor = zeros( size( b ) );
for t=s:s+TimeToEscapeBounds
    for ii=1:length(X_list)
        b( t ) = b( t ) + A( ii+1 , b_i ) * eval(char(X_list( ii )));
        add_factor( t ) = add_factor( t ) + A( ii+1 , add_factors_i ) * eval(char(X_list( ii )));
    end
end

c = ( p.b_limit-b )';
c = [ c ; -y ];
ceq = ( ( b(s:s+TimeToEscapeBounds)-p.b_limit )'.* y );

end