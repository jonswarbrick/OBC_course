%% Solves SOE model with a borrowing constraint using dynareOBC perfect foresight solver
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

%% With OBC
dynareOBC soe_borrowing_constraint shockscale=-2 timetoescapebounds=100

VarNames = char('c','h','b','mu');
ShockNames = char('epsz');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/irfs.mat','irfs','IRFoffset')

clear;

%% Without OBC
dynareOBC soe_borrowing_constraint shockscale=-2 bypass

VarNames = char('c','h','b','mu');
ShockNames = char('epsz');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/irfs_unbounded.mat','irfs','IRFoffset')

clear;

%% Cleanup
dynare_cleanup(  );


%% Plot comparison
close all;
irf_plots;
