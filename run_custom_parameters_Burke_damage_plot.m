close all
clear all

%% load

% global abatement_cost_on
% abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 
% 
% global Burke_damage_on;
% Burke_damage_on = 1;  %1 means Burke damage through TFP and capital depreciation pathways, 0 means default DICE, 2 means no damage
% 
% global optimize_on;
% optimize_on = 0;  %1 means optimizing, 0 means its reading in a 1.5
% control rate, 2 means its reading a 2.0 control rate
% 
% global optimize_for_total_utility;
% optimize_for_total_utility = 1; %1 means optimize for total utility, 0 means optimize for per-capita utility

%global utility_function;
%utility_function = 1; %1 means utility function in code and 0 means utility function in text

% ----------

load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard_a1_B1_o1_u0_uf2.mat',...
  'overall_plot_series',...
  'NPV_utility',...
  'code_string')

overall_plot_series_a1_B1_o1_u0_uf2 = overall_plot_series;
NPV_utility_a1_B1_o1_u0_uf2 = NPV_utility;

% ----------

load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard_a1_B1_o0_u1_uf1.mat',...
  'overall_plot_series',...
  'NPV_utility',...
  'code_string')

overall_plot_series_a1_B1_o0_u1_uf1 = overall_plot_series;
NPV_utility_a1_B1_o0_u1_uf1 = NPV_utility;

% ----------

load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard_a1_B1_o2_u1_uf1.mat',...
  'overall_plot_series',...
  'NPV_utility',...
  'code_string')

overall_plot_series_a1_B1_o2_u1_uf1 = overall_plot_series;
NPV_utility_a1_B1_o2_u1_uf1 = NPV_utility;

%experiments are:

%1) Burke damage, optimizing, free abatement
%2) Burke damage, optimizing, with abatement cost
%3) Burke damage, control emissions as if abatement is free but actually
%pay abatement costs
%4) DICE damage, optimizing, with abatement cost

overall_plot_series_all_exps = cat(3,overall_plot_series_a1_B1_o2_u1_uf1,...
                                     overall_plot_series_a1_B1_o0_u1_uf1,...
                                     overall_plot_series_a1_B1_o1_u0_uf2);
    
%% define

titles = {'total utility',...
          'per-capita utility',...
          'gross output',...
          'per-capita consumption',...
          'total capital',...
          'CO_2 intensity of output',...
          'population',...
          'discount factor',...
          'net output (Y_G)',...
          'industrial CO_2 emissions',...
          'ATM CO_2',...
          'global temperature',...
          'warming damage',...
          'control rate',...
          'carbon tax',...
          'mitigation cost'};
      
ylabels = {'utils',...
          'utils / person',...
          '$ trillion 2005 USD',...
          '$ thousand 2005 USD / person',...
          '$ trillion 2005 USD',...
          'MTC/$1000',...
          'million people',...
          'fraction',...
          '$ trillion 2005 USD',...
          'Gt CO_2 / year',...
          'Gt CO_2',...
          '\Delta T above PI',...
          'fraction GWP lost',...
          'fraction',...
          '$1000 / Ton CO_2',...        
          'fraction GWP lost'};

%% plot 

exp_lnstyles = {'-k','-.r','--b'};

time = 2010:5:2330;
      
row_num = 4;
col_num = 4;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',10);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
    end_year_plot = 2125;

for plot_i = 1:length(titles)
    
    subplot(row_num,col_num,plot_i)
    hold on
    grid on
    
        %plot all
        
        for exp_i = 1:3

            plot(time,squeeze(overall_plot_series_all_exps(plot_i,:,exp_i)),exp_lnstyles{exp_i},'LineWidth',1.5)
            
        end
              
        xlim([2015 end_year_plot])

        xlabel('Year')
        ylabel(ylabels{plot_i})
        title(titles{plot_i})
        
end  

%% make a tighter plot
% 
% titles = {'total utility',...               1
%           'per-capita utility',...          2 
%           'gross output',...           3
%           'per-capita consumption',...      4
%           'total capital',...               5
%           'CO_2 intensity of output',...    6
%           'population',...                  7
%           'discount factor',...             8
%           'net output (Y_G)',...            9
%           'industrial CO_2 emissions',...   10
%           'ATM CO_2',...                    11
%           'global Temperature',...          12
%           'damage from CC',...              13
%           'control rate',...                14
%           'carbon tax',...                  15
%           'abatement cost'};                16

overall_plot_series_all_exps(13,:,:) = 1-overall_plot_series_all_exps(13,:,:);

%discount everything

overall_plot_series_total_discounted = NaN(size(overall_plot_series_all_exps));

    for var_i = 1:size(overall_plot_series_all_exps,1)
        
        overall_plot_series_total_discounted(var_i,:,:) = overall_plot_series_all_exps(var_i,:,:).*overall_plot_series_all_exps(8,:,:);

    end
    
overall_plot_series_total_disc_diff = NaN(size(overall_plot_series_all_exps));    
    
%diff everything w.r.t caase [b]    
    for exp_i = 1:3

        overall_plot_series_total_disc_diff(:,:,exp_i) = (squeeze(overall_plot_series_total_discounted(:,:,exp_i)) - squeeze(overall_plot_series_total_discounted(:,:,1)));

    end    

%% plot    
    
sub_var_inds = [10 12 16 13];

text_font_size = 8;

time_horizon = 2100;
good_time_ind = find(time == time_horizon);

row_num = 3;
col_num = 2;

    FigHandle = figure('Position', [100, 100, 1150, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',8);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    
    end_year_plot = 2110;
    
    abatement_text_Xs = [2070 2070 2070];
    abatement_text_Ys = [0.16 0.13 0.1];

    damage_text_Xs = [2020 2020 2020];
    damage_text_Ys = [0.16 0.13 0.1];
    
    text_Xs = vertcat(abatement_text_Xs,damage_text_Xs);
    text_Ys = vertcat(abatement_text_Ys,damage_text_Ys);
    
for plot_i = 1:length(sub_var_inds)
    
    subplot(row_num,col_num,plot_i)
    hold on
    
       for exp_i = 1:3

            plot(time,squeeze(overall_plot_series_all_exps(sub_var_inds(plot_i),:,exp_i)),exp_lnstyles{exp_i},'LineWidth',2)
            
       end
       
        xlim([2015 end_year_plot])
        
        if plot_i == 4; ylim([0 0.2]); end        
        if plot_i == 3; ylim([0 0.2]); end
        if plot_i == 2; ylim([0.5 2.5]); end
        if plot_i == 1; ylim([0 65]); end

        if plot_i == 2;
            
            plot([time(1) end_year_plot],[1.5 1.5],'-k','LineWidth',0.5)
            plot([time(1) end_year_plot],[2 2],'-k','LineWidth',0.5)
            
        end
        
        if plot_i == 1;
            
            legend('Run [1]: Stabilize at 2C',...
                   'Run [2]: Stabilize at 1.5C',...
                   'Run [3]: Optimize NPV total utility',...
                   'Location','NorthEast'); 
        end
        
        if plot_i == 3 || plot_i == 4;
                        
            %just the weighed average of the frac lost (but this does not
            %account for the GWP growth
                NPV_case_1 = nansum(overall_plot_series_total_discounted(sub_var_inds(plot_i),1:good_time_ind,1))./...
                             nansum(overall_plot_series_all_exps(8,1:good_time_ind,1));

                NPV_case_2 = nansum(overall_plot_series_total_discounted(sub_var_inds(plot_i),1:good_time_ind,2))./...
                             nansum(overall_plot_series_all_exps(8,1:good_time_ind,2));

                NPV_case_3 = nansum(overall_plot_series_total_discounted(sub_var_inds(plot_i),1:good_time_ind,3))./....
                             nansum(overall_plot_series_all_exps(8,1:good_time_ind,3));

            
            text(text_Xs(plot_i-2,1),text_Ys(plot_i-2,1),strcat('NPV =',num2str(round(100*NPV_case_1,1)),'%'),'Color','k','FontSize',text_font_size)
            text(text_Xs(plot_i-2,2),text_Ys(plot_i-2,2),strcat('NPV =',num2str(round(100*NPV_case_2,1)),'%'),'Color','r','FontSize',text_font_size)
            text(text_Xs(plot_i-2,3),text_Ys(plot_i-2,3),strcat('NPV =',num2str(round(100*NPV_case_3,1)),'%'),'Color','b','FontSize',text_font_size)
            
                        %total discounted mitigation cost (in % of GWP)
                %take discounted GWP multiplied by abatement or damage
                %fraction. sum that
                        
        end
        
        plot([2100 2100],[0 65],'-k','LineWidth',0.5)
        
        xlabel('Year')
        ylabel(ylabels{sub_var_inds(plot_i)})
        title(titles{sub_var_inds(plot_i)})
        
end


    subplot(row_num,col_num,plot_i+1)
    hold on
    
       for exp_i = 1:3

            plot(time,squeeze(overall_plot_series_total_disc_diff(1,:,exp_i)),exp_lnstyles{exp_i},'LineWidth',1.5)
            
       end

        xlim([2015 end_year_plot])
        
            %output info to plot
            
            lost_utility_NPV_case_1 = nansum(overall_plot_series_total_disc_diff(1,1:good_time_ind,1));
            lost_utility_NPV_case_2 = nansum(overall_plot_series_total_disc_diff(1,1:good_time_ind,2));            
            lost_utility_NPV_case_3 = nansum(overall_plot_series_total_disc_diff(1,1:good_time_ind,3));
            
            %text(2060,15,strcat('sum =',num2str(round(lost_utility_NPV_case_1,0))),'Color','k','FontSize',text_font_size)
            text(2020,-75,strcat('sum =',num2str(round(lost_utility_NPV_case_2,0))),'Color','r','FontSize',text_font_size)
            text(2055,30,strcat('sum =',num2str(round(lost_utility_NPV_case_3,0))),'Color','b','FontSize',text_font_size)
            
            plot([2100 2100],[-275 55],'-k','LineWidth',0.5)
            
            ylim([-275 55])
        
        xlabel('Year')
        ylabel('utils')
        title('diff in utility w.r.t. to run [1]')
    
        
    subplot(row_num,col_num,plot_i+2)
    hold on

       for exp_i = 1:3

            plot(time,squeeze(overall_plot_series_total_disc_diff(9,:,exp_i)),exp_lnstyles{exp_i},'LineWidth',1.5)
            
       end

        xlim([2015 end_year_plot])
        
            %output info to plot
            
            lost_output_NPV_case_1 = nansum(overall_plot_series_total_disc_diff(9,1:good_time_ind,1));
            lost_output_NPV_case_2 = nansum(overall_plot_series_total_disc_diff(9,1:good_time_ind,2));            
            lost_output_NPV_case_3 = nansum(overall_plot_series_total_disc_diff(9,1:good_time_ind,3));
            
            %text(2060,0.5,strcat('sum =',num2str(round(lost_output_NPV_case_1,0))),'Color','k','FontSize',text_font_size)
            text(2020,-2.5,strcat('sum =',num2str(round(lost_output_NPV_case_2,0))),'Color','r','FontSize',text_font_size)
            text(2041,0.6,strcat('sum =',num2str(round(lost_output_NPV_case_3,0))),'Color','b','FontSize',text_font_size)
            
            plot([2100 2100],[-275 55],'-k','LineWidth',0.5)
            
            ylim([-6 3])
        
        xlabel('Year')
        ylabel('$ trillion 2005 USD')
        title('diff in GWP w.r.t. to run [1]')
        
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/test_burke_fig','-dpdf')        

    
    
    
    
    
    
    
    