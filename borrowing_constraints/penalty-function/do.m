%% Solves SOE model with a penalty function approxamation to a borrowing constraint
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

%% With OBC
dynareOBC soe_obc shockscale=-2

[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/irfs.mat','irfs','IRFoffset')

clear;

%% Cleanup
dynare_cleanup(  );




