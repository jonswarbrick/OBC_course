%% Plots IRFS simulations for multiple models
% For the course "Occasionally Binding Constraints in Macro"
% Jonathan Swarbrick, 2021
clear;

%% Some options

opt.periods_to_plot = 40;
opt.no_rows_sub_plots = 1;
opt.no_cols_sub_plots = 3;
opt.plot_size = [100 100 900 200];

opt.results_files = {
    'results/irfs_unbounded'
    'results/irfs'
};

opt.model_names = {
    'Perfect Foresight - always binding'
    'Perfect Foresight - OBC'
};

opt.variable_names = char( ...
    'c' , 'h' , 'b' ...
);

opt.variable_labels = char( ...
    'Consumption','Hours','Bonds / Cons' ...
);

opt.shock_names = char( ...
    'epsz' ...
);

opt.shock_labels = char( ...
    'Technology shock' ...
);

opt.xlabel = 'Quarters';
opt.ylabel = 'Level deviation from SS';
opt.nolegend = 1;

IRF_plotter( opt );
