%------------------------------------------------------------%
%	Switching version of									 %
%	"Global Dynamics at the Zero Lower Bound" by			 %
%	William T. Gavin, Benjamin D. Keen, Alexander W. Richter %
%	And Nathaniel A. Throckmorton							 %
%   (Model 2: Baseline with Capital)						 %
%------------------------------------------------------------%

endogenous W "Real wages"
N "Labor Hours"
C "Consumption"
R "Nominal interest rate"
BETA "Discount factor"
PAI	"Inflation rate"
Z "TFP process"
Y "Output"
PSI	"Real Marginal Cost"
RK "Rental Rate"
I "Investment"
K "Capital"
RTAYLOR

exogenous  EZ "TFP shock" EBETA	"Discount shock"

parameters	chi, eta, sigma, betass, delta, alpha, varphi, theta, paistar, phi_pai, phi_y
zbar, rhoz, rhobeta, sigmabeta, sigmaz,nu, g_over_y, g, rbar
gam	kappa1 kappa2


parameters	k_over_y n_over_y i_over_y c_over_y

@#if exogenous_switching
	zlb_tp_1_2, zlb_tp_2_1
@#end

parameters(zlb,2) zlb_flag	% indicator 

model
	# paibar = paistar;
	# omicron1 = 1/kappa1-1;
	# omicron2 = 1/kappa2-1;
	% endogenous switching imposes some disciplin:
	% you have to say how exactly you move from one regime
	% to the other and vice-versa
	%-----------------------------------------------------
	@#if ~exogenous_switching
		! zlb_tp_1_2 = 1/(1+omicron1*exp(-gam*(BETA-betass)));% 0.1;
		! zlb_tp_2_1 = 1/(1+omicron2*exp(gam*(BETA-betass)));% 0.2;
	@#end
	
	W = chi*N^eta*C^sigma;

	1 = R*BETA{+1}*(C/C{+1})^sigma/PAI{+1};

	1 = betass*(C/C{+1})^sigma*(RK{+1}+1-delta);

	Y = Z*K{-1}^alpha*N^(1-alpha);

	K = (1-delta)*K{-1}+I;

	PSI =W^(1-alpha)*RK^alpha*(1-alpha)^(1-alpha)*alpha^alpha/Z;

	K{-1}/N = alpha/(1-alpha)*W/RK;

	varphi*(PAI/paibar-1)*PAI/paibar = (1-theta)+theta*PSI+varphi*BETA{+1}*(C/C{+1})^sigma*(PAI{+1}/paibar-1)*PAI{+1}/paibar*Y{+1}/Y;

	C + I + g = (1-0.5*varphi*(PAI/paibar-1)^2)*Y;

	RTAYLOR = steady_state(R)*(PAI/steady_state(PAI))^phi_pai*(Y/steady_state(Y))^phi_y;

	R = zlb_flag*1+(1-zlb_flag)*RTAYLOR;  % <---	R = max(1,RTAYLOR);

	Z = zbar*(Z{-1}/zbar)^rhoz*exp(sigmaz*EZ);

	BETA = betass*(BETA{-1}/betass)^rhobeta*exp(sigmabeta*EBETA);

steady_state_model
	PSI = (theta-1)/theta;
	Z = zbar;
	BETA =  betass;
	RK = 1/BETA-1+delta;
	rbar=paistar/BETA;
	RTAYLOR = rbar;
	R = zlb_flag*1+(1-zlb_flag)*RTAYLOR;	
	PAI = BETA*R;
	W = (PSI*Z/(RK^alpha*(1-alpha)^(1-alpha)*alpha^alpha))^(1/(1-alpha));
	i_over_y = delta/Z*(alpha/(1-alpha)*W/RK)^(1-alpha);%xx_ssmdef_1
	c_over_y = (1-g_over_y-i_over_y);
	n_over_y = ((1-alpha)/alpha*RK/W)^alpha/Z;
	Y = (W/(chi*n_over_y^eta*c_over_y^sigma))^(1/(eta+sigma)); 
	N = n_over_y*Y;
	I = i_over_y*Y; 
	K = alpha/(1-alpha)*W/RK*N;
	C = c_over_y*Y;
	g = g_over_y*Y;

parameterization
	chi, 4.507;
	eta, 1/3;
	sigma, 1;
	betass,0.99;
	delta, 0.025;
	alpha, 0.33;
	varphi, 58.25;
	theta, 6;
	nu, 2.5;
	g_over_y,0.2;
	paistar, 1.005;
	phi_pai, 1.5;
	phi_y, 0.125;
	zbar, 1;
	rhoz, 0.8;
	rhobeta,0.8;
	sigmabeta,0.0025;
	sigmaz, 0.012;
	zlb_flag(zlb,1),0;
	zlb_flag(zlb,2),1;
	@#if exogenous_switching
		zlb_tp_1_2, 0.1;
		zlb_tp_2_1,	0.2;
	@#end
	gam, 10;
	kappa1, 0.1;
	kappa2, 0.2;
