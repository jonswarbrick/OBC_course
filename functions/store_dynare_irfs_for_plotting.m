function [irfs,IRFoffset] = store_dynare_irfs_for_plotting( oo_, M_ , options_ )
%STORE_DYNARE_IRFS_FOR_PLOTTING stores IRFs of interest and saves to mat
%file in format for plotting file IRF_plotter.m
%Returns structs irf. and IRFoffset.:
%    - irf. contains the irf in level deviation
%    - IRFoffset contains the steady state level to add to the irf for levels
%The inputs are:
%    - oo_. dynare solutions
%    - M_. dynare model info
%    - options_. dynare options

num_shocks = M_.exo_nbr;
num_vars = M_.endo_nbr;

for j=1:num_shocks
    for i=1:num_vars
        VarName = strtrim( M_.endo_names(i,:) );
        ShockName = strtrim( M_.exo_names(j,:) );
        if options_.order==1 % order 1 - steady state
            IRFoffset.(VarName).(ShockName) = oo_.dr.ys(strmatch(VarName,M_.endo_names,'exact'));
        elseif options_.order==2 % order 2 - steady state + correction term
            IRFoffset.(VarName).(ShockName) = oo_.dr.ys(strmatch(VarName,M_.endo_names,'exact')) + oo_.dr.ghs2(oo_.dr.order_var==strmatch(VarName,M_.endo_names,'exact'));
        elseif options_.order==3 % order 3 - steady state + correction term
            IRFoffset.(VarName).(ShockName) = oo_.dr.ys(strmatch(VarName,M_.endo_names,'exact')) + oo_.dr.g_0x(oo_.dr.order_var==strmatch(VarName,M_.endo_names,'exact'));
        end
        irfs.(VarName).(ShockName) = oo_.irfs.([VarName,'_',ShockName]);
    end
end

end
