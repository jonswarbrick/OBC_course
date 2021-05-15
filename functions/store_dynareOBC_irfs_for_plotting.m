function [irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames )
%STORE_DYNAREOBC_IRFS_FOR_PLOTTING stores IRFs of interest and saves to mat
%file in format for plotting file IRF_plotter.m
%Returns structs irf. and IRFoffset.:
%    - irf. contains the irf in level deviation
%    - IRFoffset contains the level to add to the irf for levels
%The inputs are:
%    - dynareOBC_. the output from dynareOBC
%    - oo_. the output from dynare
%    - VarNames a [v x 1] char array with the v variable names of interest
%    - ShockNames an [s x 1] char array with the s shock names of interest
% For the course "Occasionally Binding Constraints in DSGE Models"
% Jonathan Swarbrick, 2019

[num_shocks ~ ] = size(ShockNames);
[num_vars ~ ] = size(VarNames);

for j=1:num_shocks
    for i=1:num_vars
        VarName = strtrim( VarNames(i,:) );
        ShockName = strtrim( ShockNames(j,:) );
        IRFoffset.(VarName).(ShockName) = dynareOBC_.IRFOffsets.([VarName,'_',ShockName]);
        irfs.(VarName).(ShockName) = oo_.irfs.([VarName,'_',ShockName]);
    end
end

end
