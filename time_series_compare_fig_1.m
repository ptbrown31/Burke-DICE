close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

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
                 
load(['/home/pbrown/Documents/MATLAB/DICE/test_Burke/scale_with_temperature_Burke_damage_aopt',code_string,'.mat'],...
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
                 
load(['/home/pbrown/Documents/MATLAB/DICE/test_Burke/scale_with_temperature_Burke_damage_aopt',code_string,'.mat'],...
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

filename = '/home/pbrown/Documents/MATLAB/DICE/test_Burke/raw_SSP_net_CO2_emissions_range_GtCO2.xlsx';

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
        
%         if var_i == 2
%                 
%             temps_in_2100_default_DICE = squeeze(overall_plot_series_all_default_DICE(sub_var_inds(var_i),good_time_ind,:));
%             temps_in_2100_Burke = squeeze(overall_plot_series_all_Burke(sub_var_inds(var_i),good_time_ind,:));
%             
%             scatter(zeros(length(temps_in_2100_default_DICE),1)+2100,temps_in_2100_default_DICE,6,...
%             'MarkerFaceColor','k',...
%             'MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',0.6)
%         
%             scatter(zeros(length(temps_in_2100_Burke),1)+2100,temps_in_2100_Burke,6,...
%             'MarkerFaceColor','r',...
%             'MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',0.6)
%                 
%         end

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
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/time_series_tight_temp_targets','-dpdf','-fillpage')



      