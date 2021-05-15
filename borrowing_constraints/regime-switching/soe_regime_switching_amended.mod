%% SOE model with borrowing constraint
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

var c1 c2 b1 mu2 z p_11 p_22;

varexo epsz;

parameters betta R delta rho chi sigma b_limit theta11_0 theta11_1 theta11_2 theta22_0 theta22_2 mean_c1 mean_c2;

betta = 0.99;
R = 1/betta;
delta = 0.01;
rho = 0.9;
chi = 0.5;
sigma = 0.01;
b_limit = -0.01;
theta11_0 = 1;
theta11_1 = 0;
theta11_2 = 0; 
theta22_0 = 1;
theta22_2 = 0;
mean_c2 = (1 + (R-1)*b_limit ) / (1+chi);
mean_c1 = 1 / (1+chi);

model;

p_11 = theta11_0 + theta11_1 * b1 + theta11_2 * z;
p_22 = theta22_0 + theta22_2 * z;

c1 = ( exp( z ) + R * b1(-1) -  b1 ) / ( 1 + chi );  % consumption in regime 1 in t-1 and t 
c2 = ( exp( z ) + R * b_limit -  b_limit ) / ( 1 + chi );  % consumption in regime 2 in t-1 and t 

1/c1 = ( p_11 * 1/c1(+1) + (1-p_11) * 1/mean_c2 ) - 2 * delta * b1;
1/c2 = ( p_22 * 1/c2(+1) + (1-p_22) * 1/mean_c1 ) - 2 * delta * b_limit + mu2;

z = rho * z(-1) + sigma * epsz;

end;

steady;
check;

shocks;
    var epsz = 1;
end;

stoch_simul( order = 1, irf = 0, periods = 0 , nodisplay , nomoments , nocorr   );