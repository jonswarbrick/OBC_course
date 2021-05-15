%% Computes IRF for BBW2016 NK model
% Uses dynareOBCs perfect-foresight solve
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021
clear; close all;
dynareOBC NK_ZLB bypass

%% Cleanup
dynare_cleanup(  );
