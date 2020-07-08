close all
clear all

%% comments

%The main results from the paper can be reporduced with the data within
%this script

%% setpath for rest of specific functions

addpath('/Users/patrickbrown/Dropbox/Matlab/new_DICE/DICE_Burke_Damage_growth')

%% load

time_desc = 2015:5:2335;

%load DICE-default

abatement_cost_on = 1;
Burke_damage_on = 0;
optimize_on = 0;
optimize_for_total_utility = 1;
utility_function = 1;

code_string = ['_a',num2str(abatement_cost_on),...
                     '_B',num2str(Burke_damage_on),...
                     '_o',num2str(optimize_on),...
                     '_u',num2str(optimize_for_total_utility),...
                     '_uf',num2str(utility_function)];
                 
load(['/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/DICE/DICE_MATLAB_files/scale_with_temperature_Burke_damage_aopt',code_string,'.mat'],...
          'opt_temp_i_1',...
          'aopt_time_series',...
          'overall_plot_series_all',...
          'R',...
          'titles',...
          'ylabels',...
          'code_string',...
          'end_year_plot',...
          'end_year_calc_1',...
          'end_year_calc_2',...
          'end_year_plot_ind',...
          'end_year_calc_ind_1',...
          'end_year_calc_ind_2',...
          'two_d_plot_trimmed_1',...
          'two_d_plot_trimmed_2',...
          'two_d_plot_trimmed_discounted_1',...
          'two_d_plot_trimmed_discounted_2',...
          'integrated_net_GDP_1',...
          'integrated_net_GDP_2')
      
      overall_plot_series_all_default_DICE = overall_plot_series_all;
      
      %load Burke damage
abatement_cost_on = 1;
Burke_damage_on = 1;
optimize_on = 0;
optimize_for_total_utility = 1;
utility_function = 1;

code_string = ['_a',num2str(abatement_cost_on),...
                     '_B',num2str(Burke_damage_on),...
                     '_o',num2str(optimize_on),...
                     '_u',num2str(optimize_for_total_utility),...
                     '_uf',num2str(utility_function)];
                 
load(['/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/DICE/DICE_MATLAB_files/scale_with_temperature_Burke_damage_aopt',code_string,'.mat'],...
          'opt_temp_i_1',...
          'aopt_time_series',...
          'overall_plot_series_all',...
          'R',...
          'titles',...
          'ylabels',...
          'code_string',...
          'end_year_plot',...
          'end_year_calc_1',...
          'end_year_calc_2',...
          'end_year_plot_ind',...
          'end_year_calc_ind_1',...
          'end_year_calc_ind_2',...
          'two_d_plot_trimmed_1',...
          'two_d_plot_trimmed_2',...
          'two_d_plot_trimmed_discounted_1',...
          'two_d_plot_trimmed_discounted_2',...
          'integrated_net_GDP_1',...
          'integrated_net_GDP_2')
      
      overall_plot_series_all_Burke = overall_plot_series_all;
      
      T = 65;
      
%% make comparison figure with SSPs

%read in raw range

filename = '/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/DICE/DICE_MATLAB_files/raw_SSP_net_CO2_emissions_range_GtCO2.xlsx';

%read

    time_CO2_emission_upper = xlsread(filename,'Sheet1','B6:B17');
    CO2_emission_upper = xlsread(filename,'Sheet1','C6:C17');
    
    time_CO2_emission_lower = xlsread(filename,'Sheet1','E6:E18');
    CO2_emission_lower = xlsread(filename,'Sheet1','F6:F18');
    
    CO2_emission_upper_time_desc = interp1(time_CO2_emission_upper,CO2_emission_upper,time_desc);
    CO2_emission_lower_time_desc = interp1(time_CO2_emission_lower,CO2_emission_lower,time_desc);

    SSP_emissions_mean = (CO2_emission_upper_time_desc + CO2_emission_lower_time_desc)/2;
    SSP_emissions_range = CO2_emission_upper_time_desc - SSP_emissions_mean;      
      
%% plot tight

titles{16} = 'Cost of climate change mitigation';
titles{13} = 'Cost of climate change (damages)';
titles{12} = 'Global Temperature';
titles{10} = 'CO_2 emissions';
titles{9} = 'Gross World Product';
titles{17} = 'Gross World Product per capita';

ylabels{16} = 'Percent of Gross World Product';
ylabels{13} = 'Percent of Gross World Product';
ylabels{12} = 'Temperature above preindustrial (C)';
ylabels{9} = 'Trillions of 2010 US$';
ylabels{17} = 'Thousands of 2010 US$ per capita';

marker_colors = {'k','r'};

good_time_ind = find(time_desc == 2100);

sub_var_inds = [10 12 16 13 9 17];

y_maxs = [100 6 0.4 0.4 1800 160];
y_mins = [-25 0 0 0 0 0];

num_scens = size(overall_plot_series_all_default_DICE,3);

row_num = 3;
col_num = 2;

end_year_plot = 2140;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',11);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    
for var_i = 1:length(sub_var_inds)
    
    subplot(row_num,col_num,var_i)
    hold on
    
       if var_i == 2;
            
            plot([time_desc(1) time_desc(end)],[1.5 1.5],'-b','LineWidth',0.5)
            plot([time_desc(1) time_desc(end)],[2.0 2.0],'-g','LineWidth',0.5)
            plot([time_desc(1) time_desc(end)],[3.0 3.0],'-m','LineWidth',0.5)
        
        end
    
    if var_i == 1
        
        plot([time_desc(1) time_desc(end)],[0 0],':k','LineWidth',0.5)
        
        %plot(time_desc,SSP_emissions_mean,'-g')
        boundedline(time_desc(2:18),SSP_emissions_mean(2:18),SSP_emissions_range(2:18),'-g','alpha','transparency',0.1)
        
        ylim([-25 100])
        
    end
    
    %plot default-DICE
        for exp_i = 1:num_scens
    
            if var_i ~= 4; plot(time_desc,squeeze(overall_plot_series_all_default_DICE(sub_var_inds(var_i),:,exp_i)),...
                    '-k','LineWidth',0.5); end
            if var_i == 4; plot(time_desc,-1*squeeze(overall_plot_series_all_default_DICE(sub_var_inds(var_i),:,exp_i)-1),...
                    '-k','LineWidth',0.5); end
        
        end
        
    %plot Burke
        for exp_i = 1:num_scens
    
            if var_i ~= 4; plot(time_desc,squeeze(overall_plot_series_all_Burke(sub_var_inds(var_i),:,exp_i)),...
                    '--r','LineWidth',0.5); end
            if var_i == 4; plot(time_desc,-1*squeeze(overall_plot_series_all_Burke(sub_var_inds(var_i),:,exp_i)-1),...
                    '--r','LineWidth',0.5); end
        
        end

        plot([2100 2100],[y_mins(var_i) y_maxs(var_i)],':k','LineWidth',0.5)
        
        xlim([2015 end_year_plot])
        ylim([y_mins(var_i) y_maxs(var_i)])
        
        if var_i == 3 || var_i == 4
            
            ylim([0 0.38])
            
        end

        xlabel('Year')
        ylabel(ylabels{sub_var_inds(var_i)})
        title(titles{sub_var_inds(var_i)})

end

%print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/time_series_tight_temp_targets','-dpdf','-bestfit')
print(gcf,'/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/DICE/DICE_MATLAB_files/time_series_tight_temp_targets','-dpdf','-fillpage')

%% do the 1.5 vs 2 comparison PDV comparison

% size(overall_plot_series_all_Burke)
% ans =
%     17    65    61
%     size(overall_plot_series_all_default_DICE)
% ans =
%     17    65    61

%define time horizon

    time_horiz = 2100;
    time_horiz_ind = find(time_desc == time_horiz);
    
    overall_plot_series_all_default_DICE_cut = overall_plot_series_all_default_DICE(:,1:time_horiz_ind,:);
    overall_plot_series_all_Burke_cut = overall_plot_series_all_Burke(:,1:time_horiz_ind,:);
    
%find closest to 1.5 and 2

    default_DICE_temps_in_2100 = squeeze(overall_plot_series_all_default_DICE_cut(12,end,:));
    Burke_temps_in_2100 = squeeze(overall_plot_series_all_Burke_cut(12,end,:));
    
    [value_15_default_DICE, index_15_default_DICE] = min(abs(default_DICE_temps_in_2100 - 1.5));
    [value_20_default_DICE, index_20_default_DICE] = min(abs(default_DICE_temps_in_2100 - 2.0));

    [value_15_Burke, index_15_Burke] = min(abs(Burke_temps_in_2100 - 1.5));
    [value_20_Burke, index_20_Burke] = min(abs(Burke_temps_in_2100 - 2.0));
    
    overall_plot_series_all_default_DICE_cut_15 = overall_plot_series_all_default_DICE_cut(:,:,index_15_default_DICE);
    overall_plot_series_all_default_DICE_cut_20 = overall_plot_series_all_default_DICE_cut(:,:,index_20_default_DICE);
    
    overall_plot_series_all_Burke_cut_15 = overall_plot_series_all_Burke_cut(:,:,index_15_Burke);
    overall_plot_series_all_Burke_cut_20 = overall_plot_series_all_Burke_cut(:,:,index_20_Burke);
    
titles{16} = 'Cost of climate change mitigation';
titles{13} = 'Cost of climate change (damages)';
titles{12} = 'Global Temperature';
titles{10} = 'CO_2 emissions';
titles{9} = 'Gross World Product';
titles{17} = 'Gross World Product per capita';

ylabels{16} = 'Percent of Gross World Product';
ylabels{13} = 'Percent of Gross World Product';
ylabels{12} = 'Temperature above preindustrial (C)';
ylabels{9} = 'Trillions of 2010 US$ per year';
ylabels{17} = 'Thousands of 2010 US$ per capita';

marker_colors = {'k','r'};

good_time_ind = find(time_desc == 2100);

sub_var_inds = [10 12 16 13 9 17];

y_maxs = [100 6 0.4 0.4 800 160];
y_mins = [-25 0 0 0 0 0];

num_scens = size(overall_plot_series_all_default_DICE,3);

row_num = 3;
col_num = 2;

end_year_plot = 2105;

    FigHandle = figure('Position', [100, 100, 2000, 1500]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',11);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    
for var_i = 1:length(sub_var_inds)
    
    subplot(row_num,col_num,var_i)
    hold on
    
       if var_i == 2;
            
            plot([time_desc(1) time_desc(end)],[1.5 1.5],'-b','LineWidth',0.5)
            plot([time_desc(1) time_desc(end)],[2.0 2.0],'-g','LineWidth',0.5)
            plot([time_desc(1) time_desc(end)],[3.0 3.0],'-m','LineWidth',0.5)
        
        end
    
    if var_i == 1
        
        plot([time_desc(1) time_desc(end)],[0 0],':k','LineWidth',0.5)
        
        %plot(time_desc,SSP_emissions_mean,'-g')
        boundedline(time_desc(2:18),SSP_emissions_mean(2:18),SSP_emissions_range(2:18),'-g','alpha','transparency',0.1)
        
        ylim([-25 100])
        
    end
    
    %plot default-DICE
        for exp_i = 1:num_scens
    
            if var_i ~= 4; plot(time_desc,squeeze(overall_plot_series_all_default_DICE(sub_var_inds(var_i),:,exp_i)),...
                    '-k','LineWidth',0.5); end
            if var_i == 4; plot(time_desc,-1*squeeze(overall_plot_series_all_default_DICE(sub_var_inds(var_i),:,exp_i)-1),...
                    '-k','LineWidth',0.5); end
        
        end
        
    %plot Burke
        for exp_i = 1:num_scens
    
            if var_i ~= 4
                    plot(time_desc,squeeze(overall_plot_series_all_Burke(sub_var_inds(var_i),:,exp_i)),...
                    '--r','LineWidth',0.5); 
            end
            if var_i == 4
                    plot(time_desc,-1*squeeze(overall_plot_series_all_Burke(sub_var_inds(var_i),:,exp_i)-1),...
                    '--r','LineWidth',0.5); 
            end
        
        end
        
     %plot temp target trajectories
     
     plot(time_desc(time_desc <= 2100),squeeze(overall_plot_series_all_default_DICE_cut_15(sub_var_inds(var_i),:)),'-m','LineWidth',5);
     plot(time_desc(time_desc <= 2100),squeeze(overall_plot_series_all_default_DICE_cut_20(sub_var_inds(var_i),:)),'-m','LineWidth',5);
     plot(time_desc(time_desc <= 2100),squeeze(overall_plot_series_all_Burke_cut_15(sub_var_inds(var_i),:)),'-g','LineWidth',5);
     plot(time_desc(time_desc <= 2100),squeeze(overall_plot_series_all_Burke_cut_20(sub_var_inds(var_i),:)),'-g','LineWidth',5);
     
        if var_i == 4
            
             plot(time_desc(time_desc <= 2100),-1*squeeze(overall_plot_series_all_default_DICE_cut_15(sub_var_inds(var_i),:)-1),'-m','LineWidth',5);
             plot(time_desc(time_desc <= 2100),-1*squeeze(overall_plot_series_all_default_DICE_cut_20(sub_var_inds(var_i),:)-1),'-m','LineWidth',5);
             plot(time_desc(time_desc <= 2100),-1*squeeze(overall_plot_series_all_Burke_cut_15(sub_var_inds(var_i),:)-1),'-g','LineWidth',5);
             plot(time_desc(time_desc <= 2100),-1*squeeze(overall_plot_series_all_Burke_cut_20(sub_var_inds(var_i),:)-1),'-g','LineWidth',5);
             
        end

        plot([2100 2100],[y_mins(var_i) y_maxs(var_i)],':k','LineWidth',0.5)
        
        xlim([2015 end_year_plot])
        ylim([y_mins(var_i) y_maxs(var_i)])
        
        if var_i == 3 || var_i == 4
            
            ylim([0 0.38])
            
        end

        xlabel('Year')
        ylabel(ylabels{sub_var_inds(var_i)})
        title(titles{sub_var_inds(var_i)})

end

%% just look at temp targets and GWP

%1.5 minus 2.0

    GWP_15_minus_20_default_DICE = overall_plot_series_all_default_DICE_cut_15(9,:) - overall_plot_series_all_default_DICE_cut_20(9,:);
    GWP_15_minus_20_Burke = overall_plot_series_all_Burke_cut_15(9,:) - overall_plot_series_all_Burke_cut_20(9,:);

    GWP_15_minus_20_default_DICE_discounted = GWP_15_minus_20_default_DICE.*R(2,time_desc <= 2100);
    GWP_15_minus_20_Burke_discounted = GWP_15_minus_20_Burke.*R(2,time_desc <= 2100);
    
%1.5 and 2.0 minus no mitigation

    GWP_15_minus_no_mit_default_DICE = overall_plot_series_all_default_DICE_cut_15(9,:) - overall_plot_series_all_default_DICE_cut(9,:,end);
    GWP_15_minus_no_mit_Burke = overall_plot_series_all_Burke_cut_15(9,:) - overall_plot_series_all_Burke_cut(9,:,end);
    
    GWP_20_minus_no_mit_default_DICE = overall_plot_series_all_default_DICE_cut_20(9,:) - overall_plot_series_all_default_DICE_cut(9,:,end);
    GWP_20_minus_no_mit_Burke = overall_plot_series_all_Burke_cut_20(9,:) - overall_plot_series_all_Burke_cut(9,:,end);

    GWP_15_minus_no_mit_default_DICE_discounted = GWP_15_minus_no_mit_default_DICE.*R(2,time_desc <= 2100);
    GWP_15_minus_no_mit_Burke_discounted = GWP_15_minus_no_mit_Burke.*R(2,time_desc <= 2100);
    
    GWP_20_minus_no_mit_default_DICE_discounted = GWP_20_minus_no_mit_default_DICE.*R(2,time_desc <= 2100);
    GWP_20_minus_no_mit_Burke_discounted = GWP_20_minus_no_mit_Burke.*R(2,time_desc <= 2100);
    
    
end_year_plot = 2105;

FigHandle = figure('Position', [100, 100, 2000, 1500]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',7);
set(0,'defaultAxesFontName', 'helvetica')
hold on

%1.5 vs 2.0
    subplot(3,2,1)
    hold on

        %plot temp target trajectories

        plot(time_desc(time_desc <= 2100),overall_plot_series_all_default_DICE_cut_15(9,:),'-k','LineWidth',3);
        plot(time_desc(time_desc <= 2100),overall_plot_series_all_default_DICE_cut_20(9,:),'-b','LineWidth',3);
        plot(time_desc(time_desc <= 2100),overall_plot_series_all_Burke_cut_15(9,:),'-r','LineWidth',3);
        plot(time_desc(time_desc <= 2100),overall_plot_series_all_Burke_cut_20(9,:),'-m','LineWidth',3);

        legend('Default-DICE, 1.5C','Default-DICE, 2.0C','Burke-DICE, 1.5C','Burke-DICE, 2.0C','Location','NorthWest')

        xlim([2015 end_year_plot])
        ylim([0 800])

        xlabel('Year')
        ylabel(ylabels{9});

    subplot(3,2,3)
    hold on

        %plot temp target trajectories

        plot([2015 end_year_plot],[0 0],':k','LineWidth',2)

        h1 = plot(time_desc(time_desc <= 2100),GWP_15_minus_20_default_DICE,'-k','LineWidth',3);
        h2 = plot(time_desc(time_desc <= 2100),GWP_15_minus_20_Burke,'-r','LineWidth',3);

        legend([h1 h2],'Default-DICE, 1.5C minus 2.0C','Burke-DICE, 1.5C minus 2.0C','Location','NorthWest')

        xlim([2015 end_year_plot])
        ylim([-25 40])

        xlabel('Year')
        ylabel(ylabels{9});

    subplot(3,2,5)
    hold on

        %plot temp target trajectories

        plot([2015 end_year_plot],[0 0],':k','LineWidth',2)

        plot(time_desc(time_desc <= 2100),GWP_15_minus_20_default_DICE_discounted,'-k','LineWidth',5);
        plot(time_desc(time_desc <= 2100),GWP_15_minus_20_Burke_discounted,'-r','LineWidth',5);

        PDV_default_DICE_15_minus_20 = 5*nansum(GWP_15_minus_20_default_DICE_discounted)
        PDV_Burke_15_minus_20 = 5*nansum(GWP_15_minus_20_Burke_discounted)

        title(strcat('PDV Default-DICE, 1.5C minus 2.0C =',num2str(round(PDV_default_DICE_15_minus_20,0)),', PDV Burke-DICE, 1.5C minus 2.0C =',num2str(round(PDV_Burke_15_minus_20,0))))

        xlim([2015 end_year_plot])
        ylim([-9 9])

        xlabel('Year')
        ylabel(ylabels{9});
        
%1.5 and 2.0 vs. no mit

    subplot(3,2,2)
    hold on

        %plot temp target trajectories

        plot(time_desc(time_desc <= 2100),overall_plot_series_all_default_DICE_cut(9,:,end),'-g','LineWidth',3);
        plot(time_desc(time_desc <= 2100),overall_plot_series_all_Burke_cut(9,:,end),'-c','LineWidth',3);

        plot(time_desc(time_desc <= 2100),overall_plot_series_all_default_DICE_cut_15(9,:),'-k','LineWidth',3);
        plot(time_desc(time_desc <= 2100),overall_plot_series_all_default_DICE_cut_20(9,:),'-b','LineWidth',3);
        plot(time_desc(time_desc <= 2100),overall_plot_series_all_Burke_cut_15(9,:),'-r','LineWidth',3);
        plot(time_desc(time_desc <= 2100),overall_plot_series_all_Burke_cut_20(9,:),'-m','LineWidth',3);

        legend('Default-DICE no mitigation',...
               'Burke-DICE no Mitigation',...
               'Default-DICE, 1.5C',...
               'Default-DICE, 2.0C',...
               'Burke-DICE, 1.5C',...
               'Burke-DICE, 2.0C','Location','NorthWest')

        xlim([2015 end_year_plot])
        ylim([0 800])

        xlabel('Year')
        ylabel(ylabels{9});

    subplot(3,2,4)
    hold on

        %plot temp target trajectories

        plot([2015 end_year_plot],[0 0],':k','LineWidth',2)

        h1 = plot(time_desc(time_desc <= 2100),GWP_15_minus_no_mit_default_DICE,'-k','LineWidth',3);
        h2 = plot(time_desc(time_desc <= 2100),GWP_15_minus_no_mit_Burke,'-r','LineWidth',3);
        h3 = plot(time_desc(time_desc <= 2100),GWP_20_minus_no_mit_default_DICE,':k','LineWidth',3);
        h4 = plot(time_desc(time_desc <= 2100),GWP_20_minus_no_mit_Burke,':r','LineWidth',3);

        legend([h1 h2 h3 h4],'Default-DICE, 1.5C minus no mitigation',...
                             'Burke-DICE, 1.5C minus no mitigation',...
                             'Default-DICE, 2.0C minus no mitigation',...
                             'Burke-DICE, 2.0C minus no mitigation',...
                             'Location','NorthWest')

        xlim([2015 end_year_plot])
        ylim([-40 55])

        xlabel('Year')
        ylabel(ylabels{9});

    subplot(3,2,6)
    hold on

        %plot temp target trajectories

        plot([2015 end_year_plot],[0 0],':k','LineWidth',2)

        plot(time_desc(time_desc <= 2100),GWP_15_minus_no_mit_default_DICE_discounted,'-k','LineWidth',5);
        plot(time_desc(time_desc <= 2100),GWP_15_minus_no_mit_Burke_discounted,'-r','LineWidth',5);
        plot(time_desc(time_desc <= 2100),GWP_20_minus_no_mit_default_DICE_discounted,':k','LineWidth',5);
        plot(time_desc(time_desc <= 2100),GWP_20_minus_no_mit_Burke_discounted,':r','LineWidth',5);
        
        
        PDV_default_DICE_15_minus_no_mit = 5*nansum(GWP_15_minus_no_mit_default_DICE_discounted)
        PDV_Burke_15_minus_no_mit = 5*nansum(GWP_15_minus_no_mit_Burke_discounted)
        PDV_default_DICE_20_minus_no_mit = 5*nansum(GWP_20_minus_no_mit_default_DICE_discounted)
        PDV_Burke_20_minus_no_mit = 5*nansum(GWP_20_minus_no_mit_Burke_discounted)
        
        
        title(strcat('PDV Default-DICE, 1.5C minus no mit =',num2str(round(PDV_default_DICE_15_minus_no_mit,0)),...
                   ', PDV Burke-DICE, 1.5C minus no mit =',num2str(round(PDV_Burke_15_minus_no_mit,0)),...
                   ', PDV Default-DICE, 2.0C minus no mit =',num2str(round(PDV_default_DICE_20_minus_no_mit,0)),...                     
                   ', PDV Burke-DICE, 2.0C minus no mit =',num2str(round(PDV_Burke_20_minus_no_mit,0))))
                 
        xlim([2015 end_year_plot])
        ylim([-9 9])

        xlabel('Year')
        ylabel(ylabels{9});

    
%% calculate mitigation costs and damage costs separately

load('/Users/patrickbrown/Dropbox/SJSU Weather Climate & Human Systems Lab/Projects/DICE/DICE_MATLAB_files/Run_DICE_no_damage_no_mitigation.mat',...
    'overall_plot_series_all_no_damage_no_mitigation',...
    'titles_no_damage_no_mitigation',...
    'ylabels_no_damage_no_mitigation')

baseline_GWP = overall_plot_series_all_no_damage_no_mitigation(9,1:time_horiz_ind);

disp('2100 SSP1 baseline GWP')
SSP2_2100_baseline_SSP = 539.332; %MESSAGE-GLOBIOM SSP2-Baseline https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=series

disp('scaling factor')
scale_GDP_ratio = SSP2_2100_baseline_SSP./baseline_GWP(end);

baseline_GWP = baseline_GWP*scale_GDP_ratio;

%mitigation costs
    %get fractional mitigation costs for both targets
        frac_mit_cost_for_20_default_DICE = overall_plot_series_all_default_DICE_cut_20(16,:);
        frac_mit_cost_for_15_default_DICE = overall_plot_series_all_default_DICE_cut_15(16,:);
        frac_mit_cost_for_20_Burke = overall_plot_series_all_Burke_cut_20(16,:);
        frac_mit_cost_for_15_Burke = overall_plot_series_all_Burke_cut_15(16,:);

    %get true mitigation costs

        mit_cost_for_20_default_DICE = frac_mit_cost_for_20_default_DICE.*baseline_GWP;
        mit_cost_for_15_default_DICE = frac_mit_cost_for_15_default_DICE.*baseline_GWP;
%         mit_cost_for_20_Burke = frac_mit_cost_for_20_Burke.*overall_plot_series_all_default_DICE(9,1:time_horiz_ind,index_20_Burke);
%         mit_cost_for_15_Burke = frac_mit_cost_for_15_Burke.*overall_plot_series_all_default_DICE(9,1:time_horiz_ind,index_15_Burke);
        mit_cost_for_20_Burke = frac_mit_cost_for_20_Burke.*baseline_GWP;
        mit_cost_for_15_Burke = frac_mit_cost_for_15_Burke.*baseline_GWP;

    %discounted true mitigation costs

        mit_cost_for_20_default_DICE_discounted = mit_cost_for_20_default_DICE.*R(3,time_desc <= 2100);
        mit_cost_for_15_default_DICE_discounted = mit_cost_for_15_default_DICE.*R(3,time_desc <= 2100);
        mit_cost_for_20_Burke_discounted = mit_cost_for_20_Burke.*R(3,time_desc <= 2100);
        mit_cost_for_15_Burke_discounted = mit_cost_for_15_Burke.*R(3,time_desc <= 2100);
        
%damages

    %get fractional damages for both targets
        frac_dam_cost_for_20_default_DICE = 1-overall_plot_series_all_default_DICE_cut_20(13,:);
        frac_dam_cost_for_15_default_DICE = 1-overall_plot_series_all_default_DICE_cut_15(13,:);
        frac_dam_cost_for_20_Burke = 1-overall_plot_series_all_Burke_cut_20(13,:);
        frac_dam_cost_for_15_Burke = 1-overall_plot_series_all_Burke_cut_15(13,:);

    %get true damage costs

        dam_cost_for_20_default_DICE = frac_dam_cost_for_20_default_DICE.*baseline_GWP;
        dam_cost_for_15_default_DICE = frac_dam_cost_for_15_default_DICE.*baseline_GWP;
        dam_cost_for_20_Burke = frac_dam_cost_for_20_Burke.*baseline_GWP;
        dam_cost_for_15_Burke = frac_dam_cost_for_15_Burke.*baseline_GWP;

    %discounted true damages costs

        dam_cost_for_20_default_DICE_discounted = dam_cost_for_20_default_DICE.*R(2,time_desc <= 2100);
        dam_cost_for_15_default_DICE_discounted = dam_cost_for_15_default_DICE.*R(2,time_desc <= 2100);
        dam_cost_for_20_Burke_discounted = dam_cost_for_20_Burke.*R(2,time_desc <= 2100);
        dam_cost_for_15_Burke_discounted = dam_cost_for_15_Burke.*R(2,time_desc <= 2100);

    
%plot
    FigHandle = figure('Position', [100, 100, 2000, 1500]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',10);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on

    %mitigation costs
        subplot(3,2,1)
        hold on

        plot(time_desc(time_desc <= 2100),baseline_GWP,'-b','LineWidth',2)

        plot(time_desc(time_desc <= 2100),baseline_GWP-mit_cost_for_20_default_DICE,'-k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),baseline_GWP-mit_cost_for_15_default_DICE,'-.k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),baseline_GWP-mit_cost_for_20_Burke,'-r','LineWidth',2)
        plot(time_desc(time_desc <= 2100),baseline_GWP-mit_cost_for_15_Burke,'-.r','LineWidth',2)
        
        legend('baseline GWP','Default DICE 2C','Default DICE 1.5C','Burke-DICE 2C','Burke-DICE 1.5C','Location','NorthWest')
        
        ylim([0 800])
        
        ylabel(ylabels{9});
        
        title('Mitigation costs')

        subplot(3,2,3)
        hold on
        
        plot(time_desc(time_desc <= 2100),mit_cost_for_20_default_DICE,'-k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),mit_cost_for_15_default_DICE,'-.k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),mit_cost_for_20_Burke,'-r','LineWidth',2)
        plot(time_desc(time_desc <= 2100),mit_cost_for_15_Burke,'-.r','LineWidth',2)
        
        ylabel(ylabels{9});

        subplot(3,2,5)
        hold on
        
        plot(time_desc(time_desc <= 2100),mit_cost_for_20_default_DICE_discounted,'-k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),mit_cost_for_15_default_DICE_discounted,'-.k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),mit_cost_for_20_Burke_discounted,'-r','LineWidth',2)
        plot(time_desc(time_desc <= 2100),mit_cost_for_15_Burke_discounted,'-.r','LineWidth',2)

        PDV_mit_cost_for_20_default_DICE = 5*nansum(mit_cost_for_20_default_DICE_discounted)
        PDV_mit_cost_for_15_default_DICE = 5*nansum(mit_cost_for_15_default_DICE_discounted)
        PDV_mit_cost_for_20_Burke = 5*nansum(mit_cost_for_20_Burke_discounted)
        PDV_mit_cost_for_15_Burke = 5*nansum(mit_cost_for_15_Burke_discounted)
        
        PDV_baseline_GWP = 5*nansum(baseline_GWP.*R(2,time_desc <= 2100));
        disp('fraction PDV of mitigation cost')
        
        PDV_mit_cost_for_15_Burke./PDV_baseline_GWP
        
        ylabel(ylabels{9});
        
        disp('my mitigation cost estimate of 1.5 w.r.t 2')
        PDV_mit_cost_for_15_Burke - PDV_mit_cost_for_20_Burke
        
        title(strcat('PDV Default-DICE 2.0C =',num2str(round(PDV_mit_cost_for_20_default_DICE,0)),...
                     ',PDV Default-DICE 1.5C =',num2str(round(PDV_mit_cost_for_15_default_DICE,0)),...
                     ',PDV Burke-DICE 2.0C =',num2str(round(PDV_mit_cost_for_20_Burke,0)),...
                     ',PDV Burke-DICE 1.5C =',num2str(round(PDV_mit_cost_for_15_Burke,0))));
    %damages
        subplot(3,2,2)
        hold on

        plot(time_desc(time_desc <= 2100),baseline_GWP,'-b','LineWidth',2)

        plot(time_desc(time_desc <= 2100),baseline_GWP-dam_cost_for_20_default_DICE,'-k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),baseline_GWP-dam_cost_for_15_default_DICE,'-.k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),baseline_GWP-dam_cost_for_20_Burke,'-r','LineWidth',2)
        plot(time_desc(time_desc <= 2100),baseline_GWP-dam_cost_for_15_Burke,'-.r','LineWidth',2)
        
        ylabel(ylabels{9});
        
        title('Damage costs')
        
        ylim([0 800])

        subplot(3,2,4)
        hold on

        plot(time_desc(time_desc <= 2100),dam_cost_for_20_default_DICE,'-k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),dam_cost_for_15_default_DICE,'-.k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),dam_cost_for_20_Burke,'-r','LineWidth',2)
        plot(time_desc(time_desc <= 2100),dam_cost_for_15_Burke,'-.r','LineWidth',2)
                
        ylabel(ylabels{9});

        subplot(3,2,6)
        hold on

        plot(time_desc(time_desc <= 2100),dam_cost_for_20_default_DICE_discounted,'-k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),dam_cost_for_15_default_DICE_discounted,'-.k','LineWidth',2)
        plot(time_desc(time_desc <= 2100),dam_cost_for_20_Burke_discounted,'-r','LineWidth',2)
        plot(time_desc(time_desc <= 2100),dam_cost_for_15_Burke_discounted,'-.r','LineWidth',2)
        
        ylabel(ylabels{9});

        PDV_dam_cost_for_20_default_DICE = 5*nansum(dam_cost_for_20_default_DICE_discounted)
        PDV_dam_cost_for_15_default_DICE = 5*nansum(dam_cost_for_15_default_DICE_discounted)
        PDV_dam_cost_for_20_Burke = 5*nansum(dam_cost_for_20_Burke_discounted)
        PDV_dam_cost_for_15_Burke = 5*nansum(dam_cost_for_15_Burke_discounted)
        
        disp('my Burke damage estimate of 1.5 w.r.t 2')
        PDV_dam_cost_for_15_Burke - PDV_dam_cost_for_20_Burke
        
        title(strcat('PDV Default-DICE 2.0C =',num2str(round(PDV_dam_cost_for_20_default_DICE,0)),...
                     ',PDV Default-DICE 1.5C =',num2str(round(PDV_dam_cost_for_15_default_DICE,0)),...
                     ',PDV Burke-DICE 2.0C =',num2str(round(PDV_dam_cost_for_20_Burke,0)),...
                     ',PDV Burke-DICE 1.5C =',num2str(round(PDV_dam_cost_for_15_Burke,0))));


      