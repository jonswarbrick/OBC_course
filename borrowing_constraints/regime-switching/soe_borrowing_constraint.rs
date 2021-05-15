%% SOE model with borrowing constraint
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

endogenous c h b z

exogenous epsz

parameters R delta rho chi siggma b_limit theta_1_2 theta_2_1 psi1 psi2

parameters(borrcon,2) conflag

model

    ! borrcon_tp_1_2 = theta_1_2/(theta_1_2+exp(psi1*(b-b_limit)));
    ! borrcon_tp_2_1 = theta_2_1/(theta_2_1+exp(psi2*z));

    c = ( exp( z ) + R * b{-1} -  b ) / ( 1 + chi );
    h = 1 - chi * c / exp(z);
    z = rho * z{-1} + siggma * epsz;
    (1-conflag) * ( 1/c - 1/c{+1} + 2 * delta * b ) + conflag * (b-b_limit) = 0;

    ? b-b_limit>=0;

steady_state_model

    z = 0;
    b = conflag * b_limit / ((1-conflag) * delta + conflag );
    c = ( 1 + (R-1) * b ) / ( 1 + chi );
    h = 1 - chi * c;

parameterization

    R, 1/0.99;
    delta, 0.01;
    rho, 0.9;
    chi, 0.5;
    siggma, 0.01;
    b_limit, -0.01;
	conflag(borrcon,1),0;
	conflag(borrcon,2),1;
	psi1, 500;
	psi2, -200;
	theta_1_2, 2;
	theta_2_1, 2;
