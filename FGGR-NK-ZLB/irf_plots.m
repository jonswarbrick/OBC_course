%% Plots IRFS simulations for multiple models
% For the Bank of Canada -- Carleton course "Occasionally Binding Constraints in Macroeconomics"
% Jonathan Swarbrick, 2021

%% Some options

opt.periods_to_plot = 40;
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 2;
opt.plot_size = [200 200 900 600];

opt.results_files = {
    'results/NK_ZLB_PF'
    'results/NK_PF'
};

opt.model_names = {
    'NK ZLB'
    'NK'
};

opt.variable_names = char( ...
    'y','pi','r','nu' ...
);

opt.variable_labels = char( ...
    'Output','Inflation','Nominal Interest Rate','Price Dispersion' ...
);

opt.shock_names = char( ...
    'epsilon_b' ...
);

opt.shock_labels = char( ...
    'Time Preference Shock' ...
);

IRF_plotter( opt );
