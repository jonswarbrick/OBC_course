%% Plots IRFS simulations for multiple models
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

%% Some options

opt.periods_to_plot = 40;
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 2;
opt.plot_size = [200 200 400 400];

opt.results_files = {
    'results/SW07_nozlb'
    'results/SW07_sol1'
};

opt.model_names = {
    'No ZLB'
    'First solution found'
};

opt.variable_names = char( ...
    'y_obs','c_obs','pi_obs','r_obs'...
);

opt.variable_labels = char( ...
    'Output','Consumption','Inflation','Nominal Interest Rate'...
);

opt.shock_names = char( ...
    'epsilon' ...
);

opt.shock_labels = char( ...
    'Shock' ...
);

IRF_plotter( opt );

opt.results_files = {
    'results/SW07_nozlb'
    'results/SW07_sol2'
};

opt.model_names = {
    'No ZLB'
    'Next solution'
};

IRF_plotter( opt );
