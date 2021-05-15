function [ ] = IRF_plotter( opt )
% Plots IRFS simulations for multiple models
% For the course "Occasionally Binding Constraints in DSGE Models"
% Jonathan Swarbrick, 2019
% 
% Inputs:
% opt.periods_to_plot
% opt.no_rows_sub_plots
% opt.no_cols_sub_plots
% opt.plot_size
% opt.results_files
% opt.model_names 
% opt.variable_names 
% opt.variable_labels 
% opt.shock_names 
% opt.shock_labels
% opt.xlabel
% opt.ylabel
% opt.nolegend : 1 if no legend required, 0 for legend (0 default)

if ~isfield(opt,'xlabel')
    opt.xlabel = 'Periods';
end
if ~isfield(opt,'ylabel')
    opt.ylabel = 'Deviation from SS';
end
if ~isfield(opt,'nolegend')
    opt.nolegend = 0;
end
if ~isfield(opt,'variable_labels')
    opt.variable_labels = opt.variable_names;
end
if ~isfield(opt,'shock_labels')
    opt.shock_labels = opt.shock_names;
end


% Data Prep
opt.size_var = size(opt.variable_names);
opt.num_var = opt.size_var(1);
opt.size_shock_names = size(opt.shock_names);
opt.number_shocks = opt.size_shock_names(1);
opt.num_datasets = length(opt.results_files);

index = 0;
for ii=1:opt.num_datasets
    clearvars -except ii opt data index
    index = index+1;
    load([cell2mat(strcat((opt.results_files(index,:)),'.mat'))]);
    for jj = 1:opt.number_shocks
    curr_shock = strtrim(opt.shock_names(jj,:));
    for kk=1:opt.num_var
        try
            if max(abs(irfs.(strtrim(opt.variable_names(kk,:))).(curr_shock)))>1e-9
                irfs.plot.(strtrim(opt.variable_names(kk,:)))(jj,:) = irfs.(strtrim(opt.variable_names(kk,:))).(curr_shock)(:,1:opt.periods_to_plot);%+IRFoffset.(strtrim(opt.variable_names(kk,:))).(curr_shock);
            else
                eval( strcat('irfs.',strtrim(opt.variable_names(kk,:)),'(',num2str(jj),',:) = zeros(1,opt.periods_to_plot);'));
                irfs.plot.(strtrim(opt.variable_names(kk,:)))(jj,:) = zeros(1,opt.periods_to_plot);
            end
        catch
            disp(['Problem storing IRF for ',strtrim(opt.variable_names(kk,:))])
            irfs.plot.(strtrim(opt.variable_names(kk,:)))(jj,:) = zeros(1,opt.periods_to_plot);
        end
    end
    end
    for jj = 1:opt.num_var
        curr_var = strtrim(opt.variable_names(jj,:));
        data(:,:,index,jj) = irfs.plot.(curr_var)(:,1:opt.periods_to_plot);
    end
end
clearvars -except opt data

%% Plots

for shocks = 1:opt.number_shocks  
    figure('Name',strtrim(opt.shock_labels(shocks,:)));
    for vars = 1:opt.num_var
        subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,vars),
        plot(data(shocks,1:opt.periods_to_plot,1,vars),'k-','LineWidth',2); hold on;
        if opt.num_datasets >1
        plot(data(shocks,1:opt.periods_to_plot,2,vars),'r-', 'LineWidth',2);
        end
        if opt.num_datasets >2
        plot(data(shocks,1:opt.periods_to_plot,3,vars),'b-','LineWidth',2);
        end
        if opt.num_datasets >3
        plot(data(shocks,1:opt.periods_to_plot,4,vars), 'g-','LineWidth',2);
        end
        if opt.num_datasets >4
        plot(data(shocks,1:opt.periods_to_plot,5,vars),'m:','LineWidth',2); 
        end
        if opt.num_datasets >5
        plot(data(shocks,1:opt.periods_to_plot,6,vars), 'g','LineStyle',':', 'LineWidth',1);
        end
        xlabel(opt.xlabel,'FontSize',9);
        ylabel(opt.ylabel,'FontSize',9);
        grid off
        titlename=strtrim(opt.variable_labels(vars,:));
        title(titlename,'FontSize',10)
        axis tight;
    end 
    if ~opt.nolegend
    legend(opt.model_names)
    end
    if isfield(opt,'plot_size')
    set(gcf,'position',opt.plot_size)
    end
end

end

