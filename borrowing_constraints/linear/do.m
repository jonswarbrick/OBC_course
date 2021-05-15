%% Computes IRF for SOE model without borrowing constraint
% Uses dynare's 1st order approximation
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

dynare soe
[ irfs , IRFoffset ] = store_dynare_irfs_for_plotting( oo_, M_ , options_ );
save('results/irfs.mat','irfs','IRFoffset')
clear;

%% Cleanup
dynare_cleanup(  );

