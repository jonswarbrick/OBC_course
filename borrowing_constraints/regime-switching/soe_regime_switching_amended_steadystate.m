function [ys,check] = soe_regime_switching_amended_steadystate(ys,exe)

global M_ 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fun_b = @(b1) ( (1-(theta11_0 + theta11_1 * b1)) * ( 1/(( 1 + ( R - 1 )*b1 ) / ( 1 + chi )) - 1/mean_c2 ) - 2 * delta * b1 );
b1_0 = 0;
b1 = fzero( fun_b , b1_0 );


p_22 = theta22_0;
c2 = ( 1 + ( R - 1 )*b_limit  ) / ( 1 + chi ); 
p_11 = theta11_0 + theta11_1 * b1;
c1 = ( 1 + ( R - 1 )*b1 ) / ( 1 + chi );
mu2 = (1 - p_22 ) * ( 1/c2  - 1/mean_c1 ) + 2 * delta * b_limit;
z = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
