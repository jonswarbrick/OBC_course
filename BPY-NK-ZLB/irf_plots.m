%% Plots IRFS simulations for multiple models
% For the course "Occasionally Binding Constraints in DSGE Models"
% Jonathan Swarbrick, 2019

%% Some options

opt.periods_to_plot = 20;
opt.no_rows_sub_plots = 3;
opt.no_cols_sub_plots = 1;
opt.plot_size = [200 200 500 600];

opt.results_files = {
    'results/BPY_nozlb'
    'results/BPY_lowomega'
};

opt.model_names = {
    'No ZLB'
    'Minimum ||q+My||_\infty solution'
};

opt.variable_names = char( ...
    'y','pi','i'...
);

opt.variable_labels = char( ...
    'Output','Inflation','Nominal Interest Rate'...
);

opt.shock_names = char( ...
    'e' ...
);

opt.shock_labels = char( ...
    'Shock' ...
);

IRF_plotter( opt );


opt.results_files = {
    'results/BPY_nozlb'
    'results/BPY_highomega'
};

opt.model_names = {
    'No ZLB'
    'Minimum ||y||_\infty solution'
};

IRF_plotter( opt );
