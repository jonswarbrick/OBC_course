%% Plots IRFS simulations for multiple models
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

%% Some options

opt.periods_to_plot = 40;
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 2;
opt.plot_size = [200 200 500 400];

opt.results_files = {
    'results/SW03_sol1'
    'results/SW03_sol2'
};

opt.model_names = {
    'First solution found'
    'Next solution'
};

opt.variable_names = char( ...
    'yobs','cobs','piobs','robs'...
);

opt.variable_labels = char( ...
    'Output','Consumption','Inflation','Nominal Interest Rate'...
);

opt.shock_names = char( ...
    'epsilon_b' ...
);

opt.shock_labels = char( ...
    'Shock' ...
);

IRF_plotter( opt );
