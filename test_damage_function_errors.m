close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

%% 4 runs

%run [1] assume Burke damage, get Burke damage
%run [2] assume Burke damage, get Nordhaus damage
%run [3] assume Nordhaus damage, get Burke damage
%run [4] assume Nordhaus damage, get Nordhaus damage

%load two optimized control rates

%load control rate (and everything else) for Burke damage
    load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard_a1_B1_o1_u1.mat',...
      'overall_plot_series',...
      'NPV_utility',...
      'code_string')

    overall_plot_series_a1_B1_o1_u1 = overall_plot_series;
    NPV_utility_a1_B1_o1_u1 = NPV_utility;

    aopt_Burke = overall_plot_series_a1_B1_o1_u1(14,:)';
    
%load control rate (and everything else) for Nordhaus damage
    load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard_a1_B0_o1_u1.mat',...
      'overall_plot_series',...
      'NPV_utility',...
      'code_string')

    overall_plot_series_a1_B0_o1_u1 = overall_plot_series;
    NPV_utility_a1_B0_o1_u1 = NPV_utility;

    aopt_Nordhaus = overall_plot_series_a1_B0_o1_u1(14,:)';
    
    
%% need to run the other two experiments

%run [2] assume Burke damage, get Nordhaus damage
%run [3] assume Nordhaus damage, get Burke damage

aopts = cat(2,aopt_Burke,aopt_Nordhaus);
actual_damage_codes = [0 1];

global optimize_on;
optimize_on = 0;  %1 means optimizing, 0 means its reading in a control rate

overall_plot_series_mistakes = [];

for aopti = 1:2
   
   if optimize_on == 0;
       
       global aopt
       aopt = aopts(:,aopti);
       
   end

%% decide on type of experiment 

T = 65;

global abatement_cost_on
abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 

global Burke_damage_on;
Burke_damage_on = actual_damage_codes(aopti);  %1 means Burke damage through TFP and capital depreciation pathways, 0 means default DICE, 2 means no damage

global optimize_for_total_utility;
optimize_for_total_utility = 1; %1 means optimize for total utility, 0 means optimize for per-capita utility


%% run main program

range_values_2nd_half = linspace(1,1.5,6);
range_values_1st_half = 1./(range_values_2nd_half(2:end));

range_values = horzcat(fliplr(range_values_1st_half),range_values_2nd_half);

%accordian the variables

%damage function
damage_dep_rate_coeffs_on_T = range_values*0.063;
damage_TFP_growth_coeffs_on_T = range_values*0.0025;

%Elasticity of marginal utility 
%of consumption
elasmu_vals = range_values*1.45;         % 3
elasmu_vals(elasmu_vals <= 1) = 1.001;   

% carbon intensity of output
Sigg0_vals = range_values*-0.0152;         % 4
Siggd_vals = range_values*-0.001;        % 5

%initial cost decline backstop 
%cost per period
pd_vals = range_values*0.025;            % 6

%abatement cost curve exponent
theta2_vals = range_values*2.6;          % 7

%asymptotic population
LA_vals = range_values*11500;            % 8

%total factor productivity
Ag0_vals = range_values*0.079;           % 9
Agd_vals = range_values*0.006;           % 10

%pure rate of time preference
prstp_vals = range_values*0.03;         % 11

%climate sensitivity
deltaT_vals = range_values*3.1;          % 12

num_vars = 12;

var_names = {'damage cap-dep-rate coeff on T',...
             'damage TFP coeff on T',...
             'elasticity of utility',...
             'C intensity, growth rate',...
             'C intensity, change in growth rate',...
             'backstop cost decline',...
             'abatement cost curve exponent',...
             'asymptotic population',...
             'TFP, growth rate',...
             'TFP, change in growth rate',...
             'pure rate of time preference',...
             'climate sensitivity'};

global damage_dep_rate_coeff_on_T
global damage_TFP_growth_coeff_on_T
global elasmu
global Sigg0
global Siggd
global pd
global theta2
global LA
% global Ag0
% global Agd
global prstp
global deltaT

middle_value_i = 6;

        %set everything to middle
        
                           damage_dep_rate_coeff_on_T = damage_dep_rate_coeffs_on_T(middle_value_i);
                           damage_TFP_growth_coeff_on_T = damage_TFP_growth_coeffs_on_T(middle_value_i);
                           elasmu = elasmu_vals(middle_value_i);
                           Sigg0 = Sigg0_vals(middle_value_i);
                           Siggd = Siggd_vals(middle_value_i);
                           pd = pd_vals(middle_value_i);      
                           theta2 = theta2_vals(middle_value_i);
                           LA = LA_vals(middle_value_i);    
%                            Ag0 = Ag0_vals(middle_value_i);
%                            Agd = Agd_vals(middle_value_i);    
                           prstp = prstp_vals(middle_value_i);
                           deltaT = deltaT_vals(middle_value_i);
                                                      


                                %run DICE
                                addpath(genpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth'))

                                DICE;

                                NPV_utility = fval;

                                %Calculate the Carbon Price
                                opt_price = pbacktime.*(aopt'./pai(1:length(aopt))).^(theta2-1); % Carbon price, $ per ton of CO2

                                overall_plot_series = cat(1,S(16,:),...
                                                    S(17,:),...
                                                    S(7,:),...
                                                    S(19,:),...
                                                    S(9,:),...
                                                    Sig,...
                                                    L,...
                                                    R,...
                                                    S(20,:),...
                                                    S(12,:),...
                                                    S(4,:),...
                                                    S(2,:),...
                                                    S(14,:),...
                                                    aopt',...
                                                    opt_price,...
                                                    S(10,:));
                                                
                                                overall_plot_series_mistakes = cat(3,overall_plot_series_mistakes,overall_plot_series);
                                                
end

%run [1] assume Burke damage, get Burke damage
%run [2] assume Burke damage, get Nordhaus damage
%run [3] assume Nordhaus damage, get Burke damage
%run [4] assume Nordhaus damage, get Nordhaus damage

overall_plot_series_4_runs = NaN(size(overall_plot_series_mistakes,1),size(overall_plot_series_mistakes,2),4);

overall_plot_series_4_runs(:,:,1) = overall_plot_series_a1_B1_o1_u1;
overall_plot_series_4_runs(:,:,2:3) = overall_plot_series_mistakes;
overall_plot_series_4_runs(:,:,4) = overall_plot_series_a1_B0_o1_u1;

    % S
    % 1 - regular capital X
    % 2 - GMST X
    % 3 - lower ocean temperature X
    % 4 - atmospheric CO2 X
    % 5 - upper ocean CO2 X
    % 6 - lower ocean CO2 X
    % 7 - gross output X
    % 8 - zero-carbon emissions captial (leave NaN)
    % 9 - total capital (same as 1)
    % 10 - abatement cost
    % 11 - epsilon (leave NaN)
    % 12 - Industrial CO2 emissions
    % 13 - CO2 Radiative Forcing
    % 14 - Climate Damage
    % 15 - 
    % 16 - total utility
    % 17 - per capita utility
    % 18 - consumption
    % 19 - consumption per capita
    % 20 - net output
    % 21 - deltaK
    % 22 - TFP growth
    % 23 - TFP
    
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
          'damage from warming fraction',...
          'control rate',...
          'carbon tax',...
          'abatement cost fraction'};
      
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
          'frac output w.r.t. no warming',...
          'fraction',...
          '$1000 / Ton CO_2',...        
          'fract output lost'};
      
%% plot e.g., discounted per-capita consumption      
    
%discount everything

overall_plot_series_total_discounted = NaN(size(overall_plot_series_4_runs));

    for var_i = 1:size(overall_plot_series_4_runs,1)
        
        overall_plot_series_total_discounted(var_i,:,:) = overall_plot_series_4_runs(var_i,:,:).*overall_plot_series_4_runs(8,:,:);

    end
    
overall_plot_series_all_disc_diff = NaN(size(overall_plot_series_4_runs));    
    
%diff everything w.r.t caase [b]
    for exp_i = 1:4

        overall_plot_series_all_disc_diff(:,:,exp_i) = (squeeze(overall_plot_series_total_discounted(:,:,exp_i)) - squeeze(overall_plot_series_total_discounted(:,:,4)));

    end  
    
%% make NPV diffs

time_desc = 2010:5:2330;

time_horizon = 2330;
good_time_ind = find(time_desc == time_horizon);

disc_diff_NPVs = NaN(size(overall_plot_series_all_disc_diff,1),size(overall_plot_series_all_disc_diff,3));

for var_i = 1:size(overall_plot_series_4_runs,1)
    for exp_i = 1:4
        
        disc_diff_NPVs(var_i,exp_i) = nansum(squeeze(overall_plot_series_all_disc_diff(var_i,1:good_time_ind,exp_i)));
        
    end
end

%% plot 

exp_lnstyles = {'-b','--k','-.r',':g'};
      
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
        
        for exp_i = 1:4

            plot(time_desc,squeeze(overall_plot_series_4_runs(plot_i,:,exp_i)),exp_lnstyles{exp_i},'LineWidth',1.5)
            
        end
        
        if plot_i == 5;
            
            legend('run [1] assume Burke, get Burke',...
                   'run [2] assume Burke, get Nordhaus',...
                   'run [3] assume Nordhaus, get Burke',...
                   'run [4] assume Nordhaus, get Nordhaus')
            
        end
              
        xlim([2015 end_year_plot])

        xlabel('Year')
        ylabel(ylabels{plot_i})
        title(titles{plot_i})
        
end  

%% plot refined

sub_var_inds = [10 12 16 4];

row_num = 3;
col_num = 2;

    FigHandle = figure('Position', [100, 100, 1150, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',7);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
    end_year_plot = 2200;
    
for plot_i = 1:length(sub_var_inds)
    
    subplot(row_num,col_num,plot_i)
    hold on
    grid on
    
       for exp_i = 1:4

            plot(time_desc,squeeze(overall_plot_series_4_runs(sub_var_inds(plot_i),:,exp_i)),exp_lnstyles{exp_i},'LineWidth',1.5)
            
       end
       
        xlim([2015 end_year_plot])
        
        if plot_i == 3; ylim([0 0.025]); end
        
        if plot_i == 2;
            
            plot([time_desc(1) end_year_plot],[1.5 1.5],'-k','LineWidth',1)
            plot([time_desc(1) end_year_plot],[2 2],'-k','LineWidth',1)
            
        end

        
%         if plot_i == 3;
%             
%             legend('Run [1]: Burke dam, opt, free abt',...
%                                'Run [2]: Burke dam, opt, DICE abt',...
%                                'Run [3]: Burke dam, act like free abt, DICE abt',...
%                                'Run [4]: DICE dam, opt, DICE abt',...
%                                'Location','NorthEast'); 
%         end
        
        xlabel('Year')
        ylabel(ylabels{sub_var_inds(plot_i)})
        title(titles{sub_var_inds(plot_i)})
        
end

text_font_size = 5;


    subplot(row_num,col_num,plot_i+1)
    hold on
    grid on
    
       for exp_i = 1:4

            plot(time_desc,squeeze(overall_plot_series_all_disc_diff(1,:,exp_i)),exp_lnstyles{exp_i},'LineWidth',1.5)
            
       end

        xlim([2015 end_year_plot])
        
            %output info to plot
            
            lost_utility_NPV_case_1 = nansum(overall_plot_series_all_disc_diff(1,1:good_time_ind,1));
            lost_utility_NPV_case_2 = nansum(overall_plot_series_all_disc_diff(1,1:good_time_ind,2));            
            lost_utility_NPV_case_3 = nansum(overall_plot_series_all_disc_diff(1,1:good_time_ind,3));
            
            text(2060,-30,strcat('sum =',num2str(round(lost_utility_NPV_case_1,0))),'Color','b','FontSize',text_font_size)
            text(2037,-7,strcat('sum =',num2str(round(lost_utility_NPV_case_2,0))),'Color','k','FontSize',text_font_size)
            text(2100,-20,strcat('sum =',num2str(round(lost_utility_NPV_case_3,0))),'Color','r','FontSize',text_font_size)
            
            %ylim([-275 25])
        
        xlabel('Year')
        ylabel('utils')
        title('discounted difference of utility w.r.t. to run [4]')
    
        
    subplot(row_num,col_num,plot_i+2)
    hold on
    grid on

       for exp_i = 1:4

            plot(time_desc,squeeze(overall_plot_series_all_disc_diff(9,:,exp_i)),exp_lnstyles{exp_i},'LineWidth',1.5)
            
       end

        xlim([2015 end_year_plot])
        
            %output info to plot

            lost_output_NPV_case_1 = nansum(overall_plot_series_all_disc_diff(9,1:good_time_ind,1));
            lost_output_NPV_case_2 = nansum(overall_plot_series_all_disc_diff(9,1:good_time_ind,2));            
            lost_output_NPV_case_3 = nansum(overall_plot_series_all_disc_diff(9,1:good_time_ind,3));
            
            text(2080,-4,strcat('sum =',num2str(round(lost_output_NPV_case_1,0))),'Color','b','FontSize',text_font_size)
            text(2140,2,strcat('sum =',num2str(round(lost_output_NPV_case_2,0))),'Color','k','FontSize',text_font_size)
            text(2120,-3,strcat('sum =',num2str(round(lost_output_NPV_case_3,0))),'Color','r','FontSize',text_font_size)
            
            %ylim([-9 1])
        
        xlabel('Year')
        ylabel('$ trillion 2005 USD')
        title('discounted difference of output w.r.t. to run [4]')    
        
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/damage_mistakes_paper','-dpdf')        
