%% SOE model with borrowing constraint
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

var c h b mu z;

varexo epsz;

parameters betta R delta rho chi sigma b_limit;

betta = 0.99;
R = 1/betta;
delta = 0.01;
rho = 0.9;
chi = 0.5;
sigma = 0.01;
b_limit = -0.01;

model;

c = ( exp( z ) + R * b(-1) -  b ) / ( 1 + chi );
h = 1 - chi * c / exp(z);
0 = min( mu , b-b_limit );
1/c = 1/c(+1) + mu - 2 * delta * b;
z = rho * z(-1) + sigma * epsz;

end;

steady;
check;

shocks;
    var epsz = 1;
end;

steady_state_model;
z = 0;
b = 0;
c = 1 / ( 1 + chi );
h = c; 
mu = 0;
end;


stoch_simul( order = 1, replic=500, irf = 80, periods = 10000 );