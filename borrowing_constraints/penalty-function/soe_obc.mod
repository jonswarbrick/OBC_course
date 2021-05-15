%% SOE model with borrowing constraint - penalty function approach
% For the course "Occasionally Binding Constraints in Macro"
% Jonathan Swarbrick, 2021

var c h b z;

varexo epsz;

parameters betta R delta rho chi sigma b_limit eta0 eta1 eta2 eta3 eta4;

betta = 0.99;
R = 1/betta;
delta = 0.01;
rho = 0.9;
chi = 0.5;
sigma = 0.01;
b_limit = -0.01;
eta0 =    1.0e+02 *0.000712397542143;
eta1 =  1.0e+02 *-0.053788161506214;
eta2 =   1.0e+02 *0.832047183276908;
eta3 =  1.0e+02 *-3.333873004382461;
eta4 =  1.0e+02 *-0.775995322015369;

model;

c = ( exp( z ) + R * b(-1) -  b ) / ( 1 + chi );
h = 1 - chi * c / exp(z);
1/c = 1/c(+1) - 2 * delta * b -( eta1 + 2*eta2*b + 3*eta3*b^2 + 4*eta4*b^3 );
z = rho * z(-1) + sigma * epsz;

end;

steady;
check;

shocks;
    var epsz = 1;
end;


stoch_simul( order = 3, irf = 80, periods = 10000 );