%% Computes IRF for FVetal2015 NK model using dynareOBC perfect foresight solver
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

dynareOBC NK_IRF ShockScale=10 nograph

VarNames = char('y','pi','r','nu');
ShockNames = char('epsilon_b');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/NK_ZLB_PF.mat','irfs','IRFoffset')

clear;

%% Without ZLB
dynareOBC NK_IRF bypass ShockScale=10 nograph

VarNames = char('y','pi','r','nu');
ShockNames = char('epsilon_b');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/NK_PF.mat','irfs','IRFoffset')
clear;


%% Cleanup
dynare_cleanup(  );

%% Plot comparison
close all;
irf_plots;
