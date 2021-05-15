%% SOE model with add-factors
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

var lead_muc b z add_factors;

varexo epsz
    news
    news_1
    news_2
    news_3
    news_4
    news_5
    news_6
    news_7
    news_8
    news_9
    news_10
    news_11
    news_12
    news_13
    news_14
    news_15
    news_16
    news_17
    news_18
    news_19
    news_20
    news_21
    news_22
    news_23
    news_24
    news_25
    news_26
    news_27
    news_28
    news_29
    news_30;

parameters betta R delta rho chi sigma b_limit;

betta = 0.99;
R = 1/betta;
delta = 0.01;
rho = 0.9;
chi = 0.5;
sigma = 0.01;
b_limit = -0.01;

model;

# c = ( exp( z ) + R * b(-1) -  b ) / ( 1 + chi );
# lead_c = ( exp( z(+1) ) + R * b -  b(+1) ) / ( 1 + chi );
# h = 1 - chi * c / exp(z);
lead_muc = 1/lead_c;
b =  exp( z ) + R * b(-1) -  ( 1 + chi )/( lead_muc - 2 * delta * b ) + add_factors;
z = rho * z(-1) + sigma * epsz;

add_factors = news + 
    news_1(-1) +
    news_2(-2) +
    news_3(-3) +
    news_4(-4) +
    news_5(-5) +
    news_6(-6) +
    news_7(-7) +
    news_8(-8) +
    news_9(-9) +
    news_10(-10) +
    news_11(-11) +
    news_12(-12) +
    news_13(-13) +
    news_14(-14) +
    news_15(-15) +
    news_16(-16) +
    news_17(-17) +
    news_18(-18) +
    news_19(-19) +
    news_20(-20) +
    news_21(-21) +
    news_22(-22) +
    news_23(-23) +
    news_24(-24) +
    news_25(-25) +
    news_26(-26) +
    news_27(-27) +
    news_28(-28) +
    news_29(-29) +
    news_30(-30);

end;


shocks;
    var epsz = 1;
    var news = 1;
    var news_1 = 1;
    var news_2 = 1;
    var news_3 = 1;
    var news_4 = 1;
    var news_5 = 1;
    var news_6 = 1;
    var news_7 = 1;
    var news_8 = 1;
    var news_9 = 1;
    var news_10 = 1;
    var news_11 = 1;
    var news_12 = 1;
    var news_13 = 1;
    var news_14 = 1;
    var news_15 = 1;
    var news_16 = 1;
    var news_17 = 1;
    var news_18 = 1;
    var news_19 = 1;
    var news_20 = 1;
    var news_21 = 1;
    var news_22 = 1;
    var news_23 = 1;
    var news_24 = 1;
    var news_25 = 1;
    var news_26 = 1;
    var news_27 = 1;
    var news_28 = 1;
    var news_29 = 1;
    var news_30 = 1;
end;

steady_state_model;
lead_muc = 1+ chi;
z = 0;
b = 0;
add_factors = 0;
end;


stoch_simul( order = 1, irf = 0, periods = 0 , nodisplay , nomoments , nocorr , nofunctions   );