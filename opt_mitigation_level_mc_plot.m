close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

%% load

time_desc = 2015:5:2335;

%load MC trials array
load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/MC_trials_Burke_vs_default_DICE_opt_mitigation.mat',...
      'overall_plot_series_all_all',...
      'titles',...
      'ylabels',...
      'index_pull',...
      'value_of_varying_variables',...
      'mc_var_names')
  
% mc_var_names = {'CO2 intensity of economy, growth rate',...
%                'CO2 intensity of economy, \Delta growth rate',...
%                'backstop cost decline rate',...
%                'mitigation cost curve exponent',...
%                'asymptotic population',...
%                'TFP, growth rate',...
%                'TFP, \Delta growth rate',...
%                '2XCO2 climate sensitivity'};  

num_mc_trials = size(overall_plot_series_all_all,4);

T = 65;

% overall_plot_series_all_all 
    % dimensions:
    % 1)variables
    % 2)time
    % 3)level of mitigation
    % 4)mc trial
    % 5)damage representation
    
% option to only retain the runs that have mitigation levels that span 1.5-3

% >> size(overall_plot_series_all_all)
% ans =
%     17    65    58   1000     2

good_time_ind = find(time_desc==2100);

for damage_repi = 1:2
    count = 2000;
    for mc_triali = 1:num_mc_trials

        temps_across_all_mitigation_levels = squeeze(overall_plot_series_all_all(12,good_time_ind,:,mc_triali,damage_repi));
        
        if min(temps_across_all_mitigation_levels) >= 1.5 || max(temps_across_all_mitigation_levels) <= 3 %these are bad
            
            overall_plot_series_all_all(:,:,:,mc_triali,damage_repi) = NaN;
            
            count = count -1;

        end
    end
    count
end
        
%% define time horizon and discount rates
tstep = 5;
discount_rates = [0 0.03 0.05];

ro = NaN(length(discount_rates),T);

for discount_rate_i = 1:length(discount_rates)
    ro(discount_rate_i,:) = discount_rates(discount_rate_i).*100.*ones(1,T);
end
%   Social time preference factor
R = zeros(length(discount_rates),T);
R(:,1) = 1;
for discount_rate_i = 1:length(discount_rates)
    for i = 2:T
        R(discount_rate_i,i) = R(discount_rate_i,i-1)/(1+ro(discount_rate_i,i-1)./100).^tstep;
    end
end

time_horzs = [2100 2335];

%% find all the different optimal temperature levels for each time horizone and discount rate

% >> size(overall_plot_series_all_all)
% ans =
%     17    65    58   1000     2

num_scens = size(overall_plot_series_all_all,3);

opt_mitigation_level_temp = NaN(num_mc_trials,length(discount_rates),length(time_horzs),2);
net_costs_15_2_3 = NaN(num_mc_trials,length(discount_rates),length(time_horzs),2,3);

for discounti = 1:length(discount_rates)
    for time_horzsi = 1:length(time_horzs)
        for damage_repi = 1:2
            for mc_triali = 1:num_mc_trials
            
                %pull out time segments of variable of interest in discount rate
                    time_start_ind = 1;
                    time_end_ind = find(time_desc == time_horzs(time_horzsi));

                    per_capita_GDP_time_seg = squeeze(overall_plot_series_all_all(17,time_start_ind:time_end_ind,:,mc_triali,damage_repi));
                    GDP_time_seg = squeeze(overall_plot_series_all_all(9,time_start_ind:time_end_ind,:,mc_triali,damage_repi));

                    R_time_seg = R(discounti,time_start_ind:time_end_ind)';

                %discount each time segment (for each mitigation level)
                    per_capita_GDP_time_seg_disc = NaN(size(per_capita_GDP_time_seg));
                    GDP_time_seg_disc = NaN(size(GDP_time_seg));

                    for exp_i = 1:num_scens

                        per_capita_GDP_time_series_this_exp = real(squeeze(per_capita_GDP_time_seg(:,exp_i)));
                        GDP_time_series_this_exp = real(squeeze(GDP_time_seg(:,exp_i)));

                        %apply discount
%                         per_capita_GDP_time_series_this_exp_disc = (per_capita_GDP_time_series_this_exp.*R_time_seg)./nansum(R_time_seg);
%                         GDP_time_series_this_exp_disc = (GDP_time_series_this_exp.*R_time_seg)./nansum(R_time_seg);
                        
                        per_capita_GDP_time_series_this_exp_disc = (per_capita_GDP_time_series_this_exp.*R_time_seg);
                        GDP_time_series_this_exp_disc = (GDP_time_series_this_exp.*R_time_seg);
                        
                        %store discounted time series
                        per_capita_GDP_time_seg_disc(:,exp_i) = per_capita_GDP_time_series_this_exp_disc;
                        GDP_time_seg_disc(:,exp_i) = GDP_time_series_this_exp_disc;

                    end

                    %integrate over time

                    npv_per_capita_GDP_time_seg_disc = nansum(per_capita_GDP_time_seg_disc,1);
                    npv_GDP_time_seg_disc = nansum(GDP_time_seg_disc,1);

                    %find mitigation level with optimim GDP per capita

                    ind_with_max_GDP_per_cap = find(npv_per_capita_GDP_time_seg_disc == max(npv_per_capita_GDP_time_seg_disc));
                    ind_with_max_GDP = find(npv_GDP_time_seg_disc == max(npv_GDP_time_seg_disc));
                    
                    %label mitigation levels with max temperature
                    
                    max_T = NaN(num_scens,1);

                    for exp_i = 1:num_scens

                        %max_T(exp_i) = max(squeeze(overall_plot_series_all_all(12,:,exp_i,mc_triali,damage_repi)));
                        max_T(exp_i) = overall_plot_series_all_all(12,good_time_ind,exp_i,mc_triali,damage_repi);

                    end

                    %put in the optimum mitigation level (labelled with max temp)
                    
                    opt_mitigation_level_temp(mc_triali,discounti,time_horzsi,damage_repi) = max_T(ind_with_max_GDP_per_cap(1));
                    
                    %place the net present value closest to given temp
                    %values into an array
                    
                   [mitigation_value_1_5, mitigation_index_1_5] = min(abs(max_T - 1.5));
                   [mitigation_value_2, mitigation_index_2] = min(abs(max_T - 2));
                   [mitigation_value_3, mitigation_index_3] = min(abs(max_T - 3));
                   
                   npv_GDP_time_seg_disc_anom = npv_GDP_time_seg_disc - npv_GDP_time_seg_disc(end);
                   npv_per_capita_GDP_time_seg_disc_anom = npv_per_capita_GDP_time_seg_disc - npv_per_capita_GDP_time_seg_disc(end);
                                      
                   net_costs_15_2_3(mc_triali,discounti,time_horzsi,damage_repi,1) = npv_GDP_time_seg_disc_anom(mitigation_index_1_5);
                   net_costs_15_2_3(mc_triali,discounti,time_horzsi,damage_repi,2) = npv_GDP_time_seg_disc_anom(mitigation_index_2);
                   net_costs_15_2_3(mc_triali,discounti,time_horzsi,damage_repi,3) = npv_GDP_time_seg_disc_anom(mitigation_index_3);
                   
%                    net_costs_15_2_3(mc_triali,discounti,time_horzsi,damage_repi,1) = npv_per_capita_GDP_time_seg_disc_anom(mitigation_index_1_5);
%                    net_costs_15_2_3(mc_triali,discounti,time_horzsi,damage_repi,2) = npv_per_capita_GDP_time_seg_disc_anom(mitigation_index_2);
%                    net_costs_15_2_3(mc_triali,discounti,time_horzsi,damage_repi,3) = npv_per_capita_GDP_time_seg_disc_anom(mitigation_index_3);
                    
            end
        end
    end
end

%% plot price histogram 

discounti = 2;
%time_horzsi = 2;

opt_temp_bins_1 = linspace(-100,50,25);
opt_temp_bins_2 = linspace(-100,50,25);
opt_temp_bins_3 = linspace(-100,700,25);
opt_temp_bins_4 = linspace(-100,700,25);

y_max = 0.6;

mitigation_level_colors = {'b','g','m'};

row_num = 2;
col_num = 2;

    FigHandle = figure('Position', [100, 100, 750, 500]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',11);
    set(0,'defaultAxesFontName', 'helvetica')
    
    lin_count = 1;
    
for time_horzsi = 1:2    
    for damage_repi = 1:2
   
            subplot(row_num,col_num,lin_count)
            hold on
            
            if lin_count == 1; opt_temp_bins = opt_temp_bins_1; end
            if lin_count == 2; opt_temp_bins = opt_temp_bins_2; end
            if lin_count == 3; opt_temp_bins = opt_temp_bins_3; end
            if lin_count == 4; opt_temp_bins = opt_temp_bins_4; end
            
            for mitig_leveli = 1:3

                net_cost_dist_to_plot = squeeze(net_costs_15_2_3(:,discounti,time_horzsi,damage_repi,mitig_leveli));
                net_cost_dist_to_plot_no_zeros = net_cost_dist_to_plot(net_cost_dist_to_plot ~= 0);

                histogram(net_cost_dist_to_plot_no_zeros,opt_temp_bins,...
                            'Normalization',...
                            'probability',...
                            'FaceColor',...
                            mitigation_level_colors{mitig_leveli},...
                            'EdgeColor',...
                            mitigation_level_colors{mitig_leveli},...
                            'EdgeAlpha',...
                            0.3,...
                            'FaceAlpha',...
                            0.3);
                        
%                 histogram(net_cost_dist_to_plot_no_zeros,...
%                             'Normalization',...
%                             'probability',...
%                             'FaceColor',...
%                             mitigation_level_colors{mitig_leveli},...
%                             'EdgeColor',...
%                             'none',...
%                             'EdgeAlpha',...
%                             0.5,...
%                             'FaceAlpha',...
%                             0.5);
                        
                        dist_median = nanmean(net_cost_dist_to_plot_no_zeros);
                        
                       plot([dist_median dist_median],[0 y_max],strcat('-',mitigation_level_colors{mitig_leveli}),'LineWidth',2)
         
            end
            
              plot([0 0],[0 y_max],'--k','LineWidth',2)
            
%             plot([1.5 1.5],[0 0.55],'--k')
%             plot([2 2],[0 0.55],'--k')
%             plot([3 3],[0 0.55],'--k')
            
             %xlabel('Net benefit of mitigation effort')
             xlabel('Gross World Product (trillion 2010 US$)')             
             ylabel('Relative Frequency')
             %title('break-even year distribution')

              xlim([min(opt_temp_bins) max(opt_temp_bins)])
              ylim([0 y_max])
            
            lin_count = lin_count + 1;
    end
end

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/net_cost_mc_plot','-dpdf','-bestfit')    


%% make a 1.5-2 histogram 

net_costs_15_minus_2 = net_costs_15_2_3(:,:,:,:,1) - net_costs_15_2_3(:,:,:,:,2);

diff_temp_bins_2100 = linspace(-35,35,21);
diff_temp_bins_2300 = -35:5:160;   %linspace(-35,132,31);

damage_colors = {'k','r'};

row_num = 1;
col_num = 2;

y_max = 0.4;

    FigHandle = figure('Position', [100, 100, 600, 200]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',6);
    set(0,'defaultAxesFontName', 'helvetica')
    
    plot_count = 1;
    
    %for discounti = 2:length(discount_rates)
    discounti = 2;
        for time_horzsi = 1:length(time_horzs)
            
            if time_horzsi == 1; diff_temp_bins = diff_temp_bins_2100; end
            if time_horzsi == 2; diff_temp_bins = diff_temp_bins_2300; end
            
            subplot(row_num,col_num,plot_count)
            hold on
                            
            for damage_repi = 1:2
    
                difference_dist_to_plot = squeeze(net_costs_15_minus_2(:,discounti,time_horzsi,damage_repi));
                
                %dont plot zeros (diddnt resolve difference because lowest
                %temp was above 2)
                
                difference_dist_to_plot_no_zeros = difference_dist_to_plot(difference_dist_to_plot ~=0);

                histogram(difference_dist_to_plot_no_zeros,diff_temp_bins,...
                            'Normalization',...
                            'probability',...
                            'FaceColor',...
                            damage_colors{damage_repi},...
                            'EdgeColor',...
                            damage_colors{damage_repi},...
                            'EdgeAlpha',...
                            0.2,...
                            'FaceAlpha',...
                            0.5);
                        
%                 histogram(difference_dist_to_plot_no_zeros,...
%                             'Normalization',...
%                             'probability',...
%                             'FaceColor',...
%                             damage_colors{damage_repi},...
%                             'EdgeColor',...
%                             'none',...
%                             'EdgeAlpha',...
%                             0.5,...
%                             'FaceAlpha',...
%                             0.5);
                        
                        
                        dist_median = nanmean(difference_dist_to_plot_no_zeros);
                        
                        damage_repi
                        time_horzsi
                        discounti
                        dist_median
                        
                        plot([dist_median dist_median],[0 y_max],strcat('-',damage_colors{damage_repi}),'LineWidth',2)
         
            end
            
            plot([0 0],[0 y_max],'--k','LineWidth',2)
            
             xlabel('Gross World Product (trillion 2010 US$)')
             ylabel('Relative Frequency')
             %title('break-even year distribution')

             %xlim([min(diff_temp_bins) max(diff_temp_bins)])
             %xlim([-5 18])
             ylim([0 y_max])
            
            plot_count = plot_count + 1;
        end
    %end

%print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/opt_mitigation_level_2_min_1_hist','-dpdf','-bestfit')
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/opt_mitigation_level_2_min_1_hist','-dpdf')

%% plot histogram

opt_temp_bins = linspace(0.5,5,25);

y_max = 0.19;

damage_colors = {'k','r'};

row_num = 3;
col_num = 2;

    FigHandle = figure('Position', [100, 100, 600, 500]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',9);
    set(0,'defaultAxesFontName', 'helvetica')
    
    plot_count = 1;
    
    for discounti = 1:length(discount_rates)
        for time_horzsi = 1:length(time_horzs)
            
            subplot(row_num,col_num,plot_count)
            hold on
                            
            for damage_repi = 1:2
    
                opt_temp_dist_to_plot = squeeze(opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi));

                histogram(opt_temp_dist_to_plot,opt_temp_bins,...
                            'Normalization',...
                            'probability',...
                            'FaceColor',...
                            damage_colors{damage_repi},...
                            'EdgeColor',...
                            damage_colors{damage_repi},...
                            'EdgeAlpha',...
                            0.3,...
                            'FaceAlpha',...
                            0.5);
                        
                        
                        dist_median = nanmean(opt_temp_dist_to_plot);
                        
                        damage_repi
                        time_horzsi
                        discounti
                        dist_median
                        
                        plot([dist_median dist_median],[0 y_max],strcat('-',damage_colors{damage_repi}),'LineWidth',2)
         
            end
            
            plot([1.5 1.5],[0 y_max],'-b','LineWidth',1.5)
            plot([2 2],[0 y_max],'-g','LineWidth',1.5)
            plot([3 3],[0 y_max],'-m','LineWidth',1.5)
            
             xlabel('optimal level of mitigation effort (C)')
             ylabel('relative frequency')
             %title('break-even year distribution')

             xlim([min(opt_temp_bins) max(opt_temp_bins)])
             ylim([0 y_max])
            
            plot_count = plot_count + 1;
        end
    end

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/opt_mitigation_level_mc_plot','-dpdf','-bestfit')


%% make a regression plot 

%get value of variables

value_of_mc_vars = NaN(length(mc_var_names),num_mc_trials);

for mc_varsi = 1:length(mc_var_names)
    for mc_triali = 1:num_mc_trials
        
        value_of_mc_vars(mc_varsi,mc_triali) = value_of_varying_variables(mc_varsi,index_pull(mc_varsi,mc_triali));
        
    end
end

%% find rankings in terms of importance 1st value mean minus 2nd value mean
    range_values_2nd_half = linspace(1,1.5,6);
    range_values_1st_half = 1./(range_values_2nd_half(2:end));

    range_values = horzcat(fliplr(range_values_1st_half),range_values_2nd_half);
    
    
%set discount rate, time horizon, and damage
discounti = 2;
time_horzsi = 2;

% linear_slope_value_variable = NaN(length(mc_var_names),2);
% correlation_value_variable = NaN(length(mc_var_names),2);

% x2_value = NaN(length(mc_var_names),2);
% x1_value = NaN(length(mc_var_names),2);

%x1_value = NaN(length(mc_var_names),2);

level_of_mitigation_difference = NaN(length(mc_var_names),2);

%n=2;

for damage_repi = 1:2
    for mc_varsi = 1:length(mc_var_names)
        
        %Predictand_Y = squeeze(opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi));
        %Predictor_X = squeeze(value_of_mc_vars(mc_varsi,:))';
        %Predictor_X = Predictor_X
        %Predictor_X = range_values(index_pull(mc_varsi,:))';
        
%         p = polyfit(Predictor_X,Predictand_Y,n);
%         
%         x2_value(mc_varsi,damage_repi) = p(1);
%         x1_value(mc_varsi,damage_repi) = p(2);
        
%         Mdl = fitlm(Predictor_X,Predictand_Y,'linear');
%         
%         linear_slope_value_variable(mc_varsi,damage_repi) = Mdl.Coefficients.Estimate(2);
%         
%         corr_mat = corrcoef(Predictor_X,Predictand_Y,'rows','pairwise');
%         
%         correlation_value_variable(mc_varsi,damage_repi) = corr_mat(2,1);

        %find the mitigation effort values assocaited with the 0.666 position

        two_thrids_mitigation_values_inds = find(squeeze(index_pull(mc_varsi,:)) == 1);
        three_halves_mitigation_values_inds = find(squeeze(index_pull(mc_varsi,:)) == 11);

        level_of_mitigation_difference(mc_varsi,damage_repi) = nanmean(squeeze(opt_mitigation_level_temp(two_thrids_mitigation_values_inds,discounti,time_horzsi,damage_repi))) - ...
                                         nanmean(squeeze(opt_mitigation_level_temp(three_halves_mitigation_values_inds,discounti,time_horzsi,damage_repi)));
    
        
    end
end

%[linear_slope_value_variable_sorted, linear_slope_value_variable_sorted_i ] = sort(abs(correlation_value_variable),1,'descend');
%[x2_value_variable_sorted, x2_value_variable_sorted_i ] = sort(abs(x2_value),1,'descend');
[diff_value_variable_sorted, diff_value_variable_sorted_i ] = sort(abs(level_of_mitigation_difference),1,'descend');


 %% plot

marker_colors = {'k','r'};

titles = {'a','b','c','d','e','f','g','h'}; 

mc_var_names = {'CO_2 intensity, growth rate',...
               'CO_2 intensity, change in growth rate',...
               'backstop cost decline rate',...
               'mitigation cost curve exponent',...
               'asymptotic population',...
               'TFP, growth rate',...
               'TFP, change in growth rate',...
               '2XCO_2 climate sensitivity'};
%            
% row_num = 2;
% col_num = 4;
% 
% for damage_repi = 1:2
% 
%     FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',7);
%     set(0,'defaultAxesFontName', 'helvetica')
%     
%     for mc_varsi = 1:length(mc_var_names)
%         
%         subplot(row_num,col_num,mc_varsi)
%         hold on
%         
%         correct_mc_var_index = diff_value_variable_sorted_i(mc_varsi,damage_repi);
%         
%         scatter(value_of_mc_vars(correct_mc_var_index,:),opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi),25,...
%             'MarkerFaceColor',marker_colors{damage_repi},...
%             'MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',0.2)
%         %scatter(range_values,opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi),6,'o','k','filled')
% 
% %         p = polyfit(value_of_mc_vars(correct_mc_var_index,:)',squeeze(opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi)),n);
% %         y = polyval(p,value_of_mc_vars(correct_mc_var_index,:));
% %         
% %         plot(value_of_mc_vars(correct_mc_var_index,:),y,'-r','LineWidth',2)
%         
% %         h=lsline;
% %         set(h(1),'color','r')
% %         set(h(1),'LineWidth',2)
%         
%         title(titles{mc_varsi})
%         %xlabel(mc_var_names{correct_mc_var_index})
%         title(strcat(titles{mc_varsi},',',mc_var_names{correct_mc_var_index}))
% 
%         ylabel('optimal level of mitigation effort (C)')
%         
%         x_min = min(value_of_mc_vars(correct_mc_var_index,:));
%         x_max = max(value_of_mc_vars(correct_mc_var_index,:));
%         
%         if x_min < 0
%             
%             x_min = 1.1*x_min;
%             x_max = 0.9*x_max;
%             
%         end
%         if x_min > 0
%             
%             x_min = 0.9*x_min;
%             x_max = 1.1*x_max;
% 
%         end     
%         
%         plot([x_min x_max],[1.5 1.5],'-b')
%         plot([x_min x_max],[2 2],'-g')
%         plot([x_min x_max],[3 3],'-m')
%         
%         ylim([0 5.5])
%         xlim([x_min x_max])
%         
%     end
% 
%     h=gcf;
%     set(h,'PaperOrientation','landscape');
% 
%     print(gcf,strcat('/home/pbrown/Documents/MATLAB/DICE/test_Burke/opt_mitigation_level_mc_plot_','damage_',num2str(damage_repi)),'-dpdf')
% end


%% plot on same plot
           
row_num = 2;
col_num = 4;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',7);
    set(0,'defaultAxesFontName', 'helvetica')
    
    for mc_varsi = 1:length(mc_var_names)        
        subplot(row_num,col_num,mc_varsi)
        hold on

        for damage_repi = 1:2
        
        correct_mc_var_index = diff_value_variable_sorted_i(mc_varsi,1);
        
        scatter(value_of_mc_vars(correct_mc_var_index,:),opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi),25,...
            'MarkerFaceColor',marker_colors{damage_repi},...
            'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',0.2)
        %scatter(range_values,opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi),6,'o','k','filled')

%         p = polyfit(value_of_mc_vars(correct_mc_var_index,:)',squeeze(opt_mitigation_level_temp(:,discounti,time_horzsi,damage_repi)),n);
%         y = polyval(p,value_of_mc_vars(correct_mc_var_index,:));
%         
%         plot(value_of_mc_vars(correct_mc_var_index,:),y,'-r','LineWidth',2)
        
%         h=lsline;
%         set(h(1),'color','r')
%         set(h(1),'LineWidth',2)
        
        title(titles{mc_varsi})
        %xlabel(mc_var_names{correct_mc_var_index})
        title(strcat(titles{mc_varsi},',',mc_var_names{correct_mc_var_index}))

        ylabel('optimal level of mitigation effort (C)')
        
        x_min = min(value_of_mc_vars(correct_mc_var_index,:));
        x_max = max(value_of_mc_vars(correct_mc_var_index,:));
        
        if x_min < 0
            
            x_min = 1.1*x_min;
            x_max = 0.9*x_max;
            
        end
        if x_min > 0
            
            x_min = 0.9*x_min;
            x_max = 1.1*x_max;

        end     
        
        plot([x_min x_max],[1.5 1.5],'-b')
        plot([x_min x_max],[2 2],'-g')
        plot([x_min x_max],[3 3],'-m')
        
        ylim([0 5.5])
        xlim([x_min x_max])
        
        end
    end

    h=gcf;
    set(h,'PaperOrientation','landscape');

    print(gcf,strcat('/home/pbrown/Documents/MATLAB/DICE/test_Burke/opt_mitigation_level_mc_plot_same','damage_',num2str(damage_repi)),'-dpdf')

%% compare to process-IAMs - calculate NPV of mitigation GDP loss to compar with SSP chart

num_scens = size(overall_plot_series_all_all,3);

time_2100_i = find(time_desc == 2100);

mean_temp_time_inds = find(time_desc >= 2081 & time_desc <= 2100);

%load baseline GWP

load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Run_DICE_no_damage_no_mitigation.mat',...
    'overall_plot_series_all_no_damage_no_mitigation',...
    'titles_no_damage_no_mitigation',...
    'ylabels_no_damage_no_mitigation')

damage_repi = 2;

GWP_no_mitigation = overall_plot_series_all_no_damage_no_mitigation(9,:);

GWP_with_mitigation = NaN(length(time_desc),num_scens,num_mc_trials);
mitigation_cost = NaN(length(time_desc),num_scens,num_mc_trials);

mitigation_cost_fraction = overall_plot_series_all_all(16,:,:,:,damage_repi);
%mitigation_cost_fraction(1,1,:) = 0;

for mc_triali = 1:num_mc_trials
    for exp_i = 1:num_scens

        GWP_with_mitigation(:,exp_i,mc_triali) = (1-squeeze(mitigation_cost_fraction(1,:,exp_i,mc_triali))).*squeeze(GWP_no_mitigation);
        
        mitigation_cost(:,exp_i,mc_triali) = GWP_with_mitigation(:,exp_i,mc_triali) - GWP_no_mitigation';

    end
end

mitigation_cost_discounted_5 = NaN(length(time_desc),num_scens,num_mc_trials);
GWP_no_mitigation_discounted_5 = R(discount_rate_i,:).*squeeze(GWP_no_mitigation);

mitigation_cost_frac_NPV = NaN(num_scens,num_mc_trials);
mean_temp_2081_2100 = NaN(num_scens,num_mc_trials);

discount_rate_i = 3; %3 is 5%

for mc_triali = 1:num_mc_trials
    for exp_i = 1:num_scens

        mitigation_cost_discounted_5(:,exp_i,mc_triali) = R(discount_rate_i,:)'.*squeeze(mitigation_cost(:,exp_i,mc_triali));
        
        if sum(isnan(mitigation_cost_discounted_5(:,exp_i,mc_triali))) < length(mitigation_cost_discounted_5(:,exp_i,mc_triali))

            mitigation_cost_frac_NPV(exp_i,mc_triali) = nansum(squeeze(-1*mitigation_cost_discounted_5(1:time_2100_i,exp_i,mc_triali)))./...
                                                                  nansum(GWP_no_mitigation_discounted_5(1:time_2100_i));

            mean_temp_2081_2100(exp_i,mc_triali) = mean(squeeze(overall_plot_series_all_all(12,mean_temp_time_inds,exp_i,mc_triali)));
            
        end

    end
end

% plot NPV mitigation cost vs. global temperature 

rcp_26_45_60_85_temps = [1.6 2.4 2.8 4.3];
rcp_26_45_60_85_temps_sig = [0.4 0.5 0.5 0.7];

rcp_26_45_60_85_temps_maxes = rcp_26_45_60_85_temps + rcp_26_45_60_85_temps_sig;
rcp_26_45_60_85_temps_mins = rcp_26_45_60_85_temps - rcp_26_45_60_85_temps_sig;

rcp_26_45_60_85_NVP_mitig_SSPall_maxes = (1/100)*[6.96 1.87 0.69 0.01];
rcp_26_45_60_85_NVP_mitig_SSPall_mins = (1/100)*[0.38 -0.01 0.01 -0.01];

rcp_26_45_60_85_NVP_mitig_SSPall_half_range = (rcp_26_45_60_85_NVP_mitig_SSPall_maxes - rcp_26_45_60_85_NVP_mitig_SSPall_mins)./2;
rcp_26_45_60_85_NVP_mitig_SSPall_central = rcp_26_45_60_85_NVP_mitig_SSPall_half_range+rcp_26_45_60_85_NVP_mitig_SSPall_mins;



% rcp_26_45_60_NVP_mitig_SSP5 = (1/100)*[2.34 1.28 0.69 0];
% rcp_26_45_60_NVP_mitig_SSP5_sig = (1/100)*[(2.91-2.34)/4 ...
%                                    (1.87-0.86)/4 ...
%                                    (0.69-0.56)/4 ...
%                                    0.1];
%                                
% rcp_26_45_60_NVP_mitig_SSP2 = (1/100)*[0.76 0.14 0.03];
% rcp_26_45_60_NVP_mitig_SSP2_sig = (1/100)*[(2.91-2.34)/4 ...
%                                    (1.87-0.86)/4 ...
%                                    (0.69-0.56)/4];
  
%     for mc_triali = 1:num_mc_trials
% 
%          plot(integrated_net_GDP_1_NPV_mitig(:,mc_triali),mean_temp_2081_2100(:,mc_triali),'-k','LineWidth',2)
%          %scatter(integrated_net_GDP_1_NPV_mitig(:,mc_triali),mean_temp_2081_2100(:,mc_triali),25,[0 0 0],'filled')
%      
%     end
rcp_colors_lines = {'-g','-b','-m','-r'};
rcp_colors = {'g','b','m','r'};

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',10);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    
    %plot all the MC trials
    
    plot(mitigation_cost_frac_NPV,mean_temp_2081_2100,'Color',[0.5 0.5 0.5],'LineWidth',0.01)

    
    for rcp_i = 1:length(rcp_26_45_60_85_temps)
        boundedline([rcp_26_45_60_85_NVP_mitig_SSPall_mins(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_maxes(rcp_i)],...
        [rcp_26_45_60_85_temps(rcp_i) rcp_26_45_60_85_temps(rcp_i)],...
        [rcp_26_45_60_85_temps_sig(rcp_i) rcp_26_45_60_85_temps_sig(rcp_i)],...
        rcp_colors_lines{rcp_i},'alpha','transparency',0.3)
    end
    
    for rcp_i = 1:length(rcp_26_45_60_85_NVP_mitig_SSPall_central)
        boundedline([rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i)],...
        [rcp_26_45_60_85_temps_mins(rcp_i) rcp_26_45_60_85_temps_maxes(rcp_i)],...
        [rcp_26_45_60_85_NVP_mitig_SSPall_half_range(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_half_range(rcp_i)],...
        rcp_colors_lines{rcp_i},'orientation', 'horiz','alpha','transparency',0.3)
    end    
    
%     boundedline(mean(integrated_net_GDP_1_NPV_mitig,2),...
%                 mean(mean_temp_2081_2100,2),...
%                 2.*std(mean_temp_2081_2100,[],2),...
%                 '-k','alpha','transparency',0.1)
%             
%     boundedline(mean(integrated_net_GDP_1_NPV_mitig,2),...
%                 mean(mean_temp_2081_2100,2),...
%                 2.*std(integrated_net_GDP_1_NPV_mitig,[],2),...
%                 '-k','orientation','horiz','alpha','transparency',0.1)

%             
         plot(nanmean(mitigation_cost_frac_NPV,2),nanmean(mean_temp_2081_2100,2),'-k','LineWidth',2)
         scatter(nanmean(mitigation_cost_frac_NPV,2),nanmean(mean_temp_2081_2100,2),25,[0 0 0],'filled')
         
%          plot(nanmedian(mitigation_cost_frac_NPV,2),nanmedian(mean_temp_2081_2100,2),'-k','LineWidth',2)
%          scatter(nanmedian(mitigation_cost_frac_NPV,2),nanmedian(mean_temp_2081_2100,2),25,[0 0 0],'filled')

     %just use error bar
        %         x = 1:10:100;
        %         y = [20 30 45 40 60 65 80 75 95 90];
        %         yneg = [1 3 5 3 5 3 6 4 3 3];
        %         ypos = [2 5 3 5 2 5 2 2 5 5];
        %         xneg = [1 3 5 3 5 3 6 4 3 3];
        %         xpos = [2 5 3 5 2 5 2 2 5 5];
        %         errorbar(x,y,yneg,ypos,xneg,xpos,'o')
        % 
        %      for rcp_i = 1:length(rcp_26_45_60_NVP_mitig_SSP5)
        %         errorbar(rcp_26_45_60_NVP_mitig_SSP2_central(rcp_i),...
        %                  rcp_26_45_60_temps(rcp_i),...
        %                  rcp_26_45_60_temps_sig(rcp_i),...
        %                  rcp_26_45_60_temps_sig(rcp_i),...
        %                  rcp_colors{rcp_i})
        %      end
        %      for rcp_i = 1:length(rcp_26_45_60_NVP_mitig_SSP5)
        %         errorbar(rcp_26_45_60_NVP_mitig_SSP2_central(rcp_i),...
        %                  rcp_26_45_60_temps(rcp_i),...
        %                  rcp_26_45_60_NVP_mitig_SSP2_half_range(rcp_i),...
        %                  rcp_26_45_60_NVP_mitig_SSP2_half_range(rcp_i),...
        %                  'horizontal',rcp_colors{rcp_i})
        %      end

     %scatter central values
     for rcp_i = 1:length(rcp_26_45_60_85_NVP_mitig_SSPall_central)
        scatter(rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i),rcp_26_45_60_85_temps(rcp_i),100,rcp_colors{rcp_i},'o','filled')
     end

         xlim([-0.003 0.075])
 
         xlabel('Net Present Value of Gross World Product loss')
         ylabel('Global temperature above preindustrial, 2081-2100 (C)')
         %title(titles{var_i})    

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/NPV_mitig_vs_avg_temp_2081_2100_MC','-dpdf','-bestfit')

%% axis switch
% 
%     FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',10);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
%     
%     %plot all the MC trials
%     
%     plot(mean_temp_2081_2100,mitigation_cost_frac_NPV,'Color',[0.5 0.5 0.5],'LineWidth',0.01)
% 
%     for rcp_i = 1:length(rcp_26_45_60_85_temps)
%         boundedline([rcp_26_45_60_85_temps(rcp_i) rcp_26_45_60_85_temps(rcp_i)],...
%         [rcp_26_45_60_85_NVP_mitig_SSPall_mins(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_maxes(rcp_i)],...
%         [rcp_26_45_60_85_temps_sig(rcp_i) rcp_26_45_60_85_temps_sig(rcp_i)],...
%         rcp_colors_lines{rcp_i},'orientation', 'horiz','alpha','transparency',0.3)
%     end
%     
%     for rcp_i = 1:length(rcp_26_45_60_85_NVP_mitig_SSPall_central)
%         boundedline([rcp_26_45_60_85_temps_mins(rcp_i) rcp_26_45_60_85_temps_maxes(rcp_i)],...
%         [rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i)],...
%         [rcp_26_45_60_85_NVP_mitig_SSPall_half_range(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_half_range(rcp_i)],...
%         rcp_colors_lines{rcp_i},'alpha','transparency',0.3)
%     end   
%     
%          plot(nanmean(mean_temp_2081_2100,2),nanmean(mitigation_cost_frac_NPV,2),'-k','LineWidth',2)
%          scatter(nanmean(mean_temp_2081_2100,2),nanmean(mitigation_cost_frac_NPV,2),25,[0 0 0],'filled')
% 
%      %scatter central values
%      for rcp_i = 1:length(rcp_26_45_60_85_NVP_mitig_SSPall_central)
%         scatter(rcp_26_45_60_85_temps(rcp_i),rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i),100,rcp_colors{rcp_i},'o','filled')
%      end
%      
%      plot([1.5 1.5],[-0.003 0.075],'-b')
%      plot([2 2],[-0.003 0.075],'-g')
%      plot([3 3],[-0.003 0.075],'-m')
% 
%          ylim([-0.003 0.075])
%  
%          ylabel('Net Present Value of Gross World Product loss')
%          xlabel('Temperature above preindustrial, 2081-2100 (C)')
%          %title(titles{var_i})    
% 
% print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/NPV_mitig_vs_avg_temp_2081_2100_MC_axis_switch','-dpdf','-bestfit')