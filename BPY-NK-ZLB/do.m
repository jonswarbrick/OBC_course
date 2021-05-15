%% Computes IRF for BPY2020 NK model
% Compares alternative solutions using dynareOBC
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../functions')

%% High omega
dynareOBC BPYModel FullHorizon omega=10000

VarNames = char('y','pi','i');
ShockNames = char('e');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/BPY_highomega.mat','irfs','IRFoffset')
clear;

%% Low omega
dynareOBC BPYModel FullHorizon omega=0.0001

VarNames = char('y','pi','i');
ShockNames = char('e');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/BPY_lowomega.mat','irfs','IRFoffset')
clear;

%% No ZLB
dynareOBC BPYModel bypass

VarNames = char('y','pi','i');
ShockNames = char('e');
[irfs,IRFoffset] = store_dynareOBC_irfs_for_plotting( dynareOBC_, oo_, VarNames , ShockNames );
save('results/BPY_nozlb.mat','irfs','IRFoffset')
clear;

%% Plot comparison
close all;
irf_plots;

%% Cleanup
dynare_cleanup(  );
