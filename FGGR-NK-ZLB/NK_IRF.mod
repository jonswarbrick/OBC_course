%% NK model (Fernández-Villaverde, Gordon, Guerrón-Quintana & Rubio-Ramírez (JEDC, 2015))
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

var r y g c w mc pi l pi_star nu aux1 aux2 a log_beta sg;

varexo epsilon_a epsilon_m epsilon_g epsilon_b;

parameters beta_STEADY A_STEADY Sg_STEADY, pi_STEADY y_STEADY varepsilon
theta phi_pi phi_y rho_a rho_r rho_b rho_g sigma_g sigma_b sigma_a sigma_m
vartheta psi;

beta_STEADY = 0.994;
Sg_STEADY = 0.2;
A_STEADY = 1;
pi_STEADY = log( 1.005 );
psi = 1;
vartheta = 1;
theta = 0.75;
varepsilon = 6;
phi_pi = 1.5;
phi_y = 0.25;
rho_r = 0;
rho_a = 0.9;
rho_g = 0.8;
rho_b = 0.8;
sigma_a = 0.0025;
sigma_m = 0.0025;
sigma_b = 0.0025;
sigma_g = 0.0025;
y_STEADY = log( (1 / (1 - Sg_STEADY)) * (((A_STEADY * ( ( (1 - theta * (1 / exp(pi_STEADY))^(1 - varepsilon) ) / (1 - theta) )^(1/(1 - varepsilon)) ) * ((varepsilon - 1) / varepsilon) * ( (1 - theta * beta_STEADY * exp(pi_STEADY) ^varepsilon)/(1 - theta * beta_STEADY * exp(pi_STEADY) ^(varepsilon-1)))) /(psi * ((1/(1 - Sg_STEADY)) * (( ( 1 - theta) / (1 - theta * exp(pi_STEADY) ^varepsilon) ) * ( ( (1 - theta * (1 / exp(pi_STEADY))^(1 - varepsilon) ) / (1 - theta) )^(1/(1 - varepsilon)) ) ^(-varepsilon))/ A_STEADY)^vartheta))^(1/(1 + vartheta))) );

model;
    #C = exp(c);
    #lead_C = exp(c(+1));
    #PI = exp( pi );
    #lead_PI = exp( pi(+1) );
    #beta = exp(log_beta);
    #lead_beta = exp(log_beta(+1));
    #lag_beta = exp(log_beta(-1));
    #AUX1 = exp(aux1);
    #lead_AUX1 = exp(aux1(+1));
    #AUX2 = exp(aux2);
    #lead_AUX2 = exp(aux2(+1));
    #PI_STAR = exp(pi_star);
    #lead_PI_STAR = exp(pi_star(+1));
    #R = exp(r);
    #lag_R = exp(r(-1));
    #A = exp(a);
    #lag_A = exp(a(-1));
    #Sg = exp(sg);
    #lag_Sg = exp(sg(-1));
    #Y = exp(y);
    #G = exp(g);
    #W = exp(w);
    #MC = exp(mc);
    #L = exp(l);
    #NU = exp(nu);
    #lag_NU = exp(nu(-1));

	Y = (A/NU) * L;
	G = Sg*Y;
    C = Y - G;
	W = psi * L^vartheta*C;
	MC = W/A;
	1 = R * lead_beta * ( C / lead_C ) / lead_PI;
	varepsilon * AUX1 = (varepsilon - 1) * AUX2;
	AUX1 = MC * (Y/C) + theta * lead_beta * lead_PI^(varepsilon) * lead_AUX1;
	AUX2 = PI_STAR * ((Y/C) + theta * lead_beta * ((lead_PI^(varepsilon-1))/lead_PI_STAR) * lead_AUX2);

    r = max(0, rho_r*r(-1) + (1 - rho_r)*( pi_STEADY - log(beta_STEADY) + phi_pi*( pi - pi_STEADY ) + phi_y*( y - y_STEADY) ) - sigma_m * epsilon_m);
    
	1 = theta * (PI^(varepsilon-1)) + (1 - theta) * PI_STAR^(1 - varepsilon);
	NU = theta * (PI^varepsilon) * lag_NU + (1 - theta) * PI_STAR^(-varepsilon);

	A = A_STEADY^(1 - rho_a) * lag_A ^(rho_a) * exp(sigma_a * epsilon_a);
	beta = beta_STEADY^(1 - rho_b) * (lag_beta^rho_b) * exp(sigma_b*epsilon_b);
	Sg = Sg_STEADY^(1 - rho_g) * lag_Sg^rho_g * exp(-sigma_g*epsilon_g);
end;

steady_state_model;
	A_ = A_STEADY;
    Sg_ = Sg_STEADY;
    beta_ = beta_STEADY;
    pi = pi_STEADY;
    PI_STAR_ = ( (1 - theta * (1 / exp(pi))^(1 - varepsilon) ) / (1 - theta) )^(1/(1 - varepsilon));
    NU_ = ( ( 1 - theta) / (1 - theta * exp(pi) ^varepsilon) ) * PI_STAR_ ^(-varepsilon);
    W_ = A_STEADY * PI_STAR_ * ((varepsilon - 1) / varepsilon) * ( (1 - theta * beta_ * exp(pi) ^varepsilon)/(1 - theta * beta_ * exp(pi) ^(varepsilon-1)));
    C_ = (W_ /(psi * ((1/(1 - Sg_)) * NU_/ A_STEADY)^vartheta))^(1/(1 + vartheta));
    Y_ = (1 / (1 - Sg_)) * C_;
    G_ = Sg_ * Y_;
    L_ = Y_ *NU_ / A_STEADY;
    MC_ = W_ / A_STEADY;
    R_ = exp(pi) / beta_;
    AUX1_ = W_ / A_STEADY * (Y_ /C_)/(1 - theta * beta_ * exp(pi) ^varepsilon);
    AUX2_ = PI_STAR_ * (Y_ /C_)/(1 - theta * beta_ * exp(pi) ^(varepsilon-1));

    c = log(C_);
    y = log(Y_);
    g = log(G_);
    w = log(W_);
    r = log(R_);
    a = log(A_);
    nu = log(NU_);
    mc = log(MC_);
    l = log(L_);
    sg = log(Sg_);
    log_beta = log(beta_);
    aux1 = log(AUX1_);
    aux2 = log(AUX2_);
    pi_star = log(PI_STAR_);
end;

shocks;
    var epsilon_a = 1;
    var epsilon_g = 1;
    var epsilon_b = 1;
    var epsilon_m = 1;
end;

steady;
check;

stoch_simul( order = 1, irf = 40, periods = 0, irf_shocks = ( epsilon_b ) );
