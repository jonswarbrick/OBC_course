%% SOE model with borrowing constraint
% For the course "Occasionally Binding Constraints in Macro"
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
    var epsz; periods 1;  values -2;
end;

steady_state_model;
z = 0;
b = 0;
c = 1 / ( 1 + chi );
h = c; 
mu = 0;
end;


simul( periods=400 );

for i=1:M_.endo_nbr
    VarName = strtrim( M_.endo_names(i,:) );
    irfs.(VarName).epsz = oo_.endo_simul(i,2:end)-oo_.steady_state(i);
    IRFoffset.(VarName).epsz = oo_.steady_state(i);
end

figure_per_shock=floor((M_.endo_nbr-1)/9)+1;
for i=1:M_.endo_nbr
    VarName = strtrim( M_.endo_names(i,:) );
    j=floor((i-1)/9)+1;
    p=i-(j-1)*9;
    figure( figure_per_shock*(1-1)+j );
    subplot(3,3,p)
    plot( oo_.endo_simul(i,1:20) );
    title( VarName );
    if p==1
        legend('epsz')
    end
end