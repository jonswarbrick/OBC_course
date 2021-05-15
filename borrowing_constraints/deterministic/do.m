%% Computes IRF for SOE model with and without borrowing constraint
% Uses dynare's built in perfect-foresight solve
% For the course "Occasionally Binding Constraints in Macro"
% Jonathan Swarbrick, 2021
clear; close all;
addpath('../../functions')

%% With OBC
dynare soe_obc
save('results/irfs.mat','irfs','IRFoffset')
clear;

%% Without OBC
dynare soe_unconstrained
save('results/irfs_unbounded.mat','irfs','IRFoffset')
clear;

%% Cleanup
dynare_cleanup(  );

%% Plot comparison
irf_compare;
