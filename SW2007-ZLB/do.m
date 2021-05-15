%% Computes IRF for SW2007 NK model
% Compares alternative solutions using dynareOBC
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../functions')

%% No bound
dynareOBC SW07IRF bypass

VarNames = char('y_obs','c_obs','pi_obs','r_obs');
ShockNames = char('epsilon');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/SW07_nozlb.mat','irfs','IRFoffset')
clear;

%% Solution 1
dynareOBC SW07IRF

VarNames = char('y_obs','c_obs','pi_obs','r_obs');
ShockNames = char('epsilon');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/SW07_sol1.mat','irfs','IRFoffset')
clear;

%% Solution 2
dynareOBC SW07IRF skipfirstsolutions=1

VarNames = char('y_obs','c_obs','pi_obs','r_obs');
ShockNames = char('epsilon');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/SW07_sol2.mat','irfs','IRFoffset')
clear;

%% Cleanup
dynare_cleanup(  );

%% Plot comparison
close all;
irf_plots;
