%% SOE model with binding borrowing constraint (regime 2)
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
1/c = 1/c(+1) + mu - delta * b;
b = b_limit;
z = rho * z(-1) + sigma * epsz;

end;

shocks;
    var epsz = 1; 
end;