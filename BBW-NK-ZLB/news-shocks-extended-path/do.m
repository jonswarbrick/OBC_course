%% Computes time-series simulation for BBW2016 NK model
% Uses dynareOBCs perfect-foresight solver
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
dynareOBC NK_ZLB timetoescapebounds=64

%% Cleanup
dynare_cleanup(  );
