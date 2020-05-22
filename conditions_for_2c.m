close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

%% decide on type of experiment 

T = 65;

global abatement_cost_on
abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 

global Burke_damage_on;
Burke_damage_on = 1;  %1 means Burke damage through TFP and capital depreciation pathways, 
                      %0 means default DICE, 
                      %2 means no damage
                      %3 means damage with a threshold for 2C so it won't go
                      %over that value
                      
global optimize_on;
optimize_on = 0;  %1 means optimizing, 0 means its reading in a control rate

%     if optimize_on == 0;
%         
         global aopt
%         aopt = ones(T,1);
%         
%     end

global optimize_for_total_utility;
optimize_for_total_utility = 1; %1 means optimize for total utility, 0 means optimize for per-capita utility

global utility_function;
utility_function = 1; %1 means utility function in code and 0 means utility function in text, 2 means utility linear with per-capita consumption

%% make a bunch of aopts to feed the model

time = 1:T;
tstep = 5; %  (5 years per period)

%make as many control rates as there are timesteps: 1 reaching 1 in each
%timestep

aopts = ones(T,T)+0.27;

for aopti = 1:T
   
    aopts(1:aopti,aopti) = linspace(0,1.27,aopti);
    
end

aopts(:,end) = zeros(T,1); %make the last one no climate policy

%revove the 1st 8 that are too fast
aopts = aopts(:,5:end);

num_scens = size(aopts,2);

%% run main program

range_values_2nd_half = linspace(1,1.5,6);
range_values_1st_half = 1./(range_values_2nd_half(2:end));

range_values = horzcat(fliplr(range_values_1st_half),range_values_2nd_half);

%accordian the variables

%damage function
damage_dep_rate_coeffs_on_T = range_values*0.013;
damage_TFP_growth_coeffs_on_T = range_values*0.0055;

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
Ag0_vals = range_values*0.076;           % 9
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
global Ag0
global Agd
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
                           Ag0 = Ag0_vals(middle_value_i);
                           Agd = Agd_vals(middle_value_i);    
                           prstp = prstp_vals(middle_value_i);
                           deltaT = deltaT_vals(middle_value_i);
                           
                           overall_plot_series_all = [];
                           
                           for aopti = 1:num_scens
                               aopti
                               
                               aopt = aopts(:,aopti);

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
                                                
                                                overall_plot_series_all = cat(3,overall_plot_series_all,overall_plot_series);
                                                
                           end
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
          'fract output w.r.t. no abatement'};
      
% Net output: 9 ($ trillion 2005 USD)
% Population: 7 (million people)
% Output/capita = $ million/person
   
per_capita_GDP_million_per_person = squeeze(overall_plot_series_all(9,:,:)./overall_plot_series_all(7,:,:));

%convert to thousand / person

per_capita_GDP_thousand_per_person = (1000).*per_capita_GDP_million_per_person;

%% get a max temp for each mitigation level

time_desc = 2015:5:2335;

max_T = NaN(num_scens,1);
   
% for exp_i = 1:num_scens
% 
%     max_T(exp_i) = max(squeeze(overall_plot_series_all(12,:,exp_i)));
% 
% end

good_time_ind = find(time_desc == 2100);
    for exp_i = 1:num_scens

        max_T(exp_i) = overall_plot_series_all(12,good_time_ind,exp_i);

    end    

%% make the discount factors

num_discount_rates = 100;
min_discount_rate = 0;
max_discount_rate = 0.08;
discount_rates = linspace(min_discount_rate,max_discount_rate,num_discount_rates);

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

% R <- (rate,time)
      
%% scan optimum temperature levels for various discount rates and time horizons 

time_horzs = 2050:5:2335;
num_time_horzs = length(time_horzs);
min_time_horz = min(time_horzs);
max_time_horz = max(time_horzs);

%empy array that will eventually be the plot
optimim_level_of_mitigation_stab_temp = NaN(length(discount_rates),length(time_horzs));

%loop through the discount rates and time horizons and calculate the best
%mitigation level for that

%per_capita_GDP_thousand_per_person <- (time,mitigation level)

for discounti = 1:num_discount_rates
    for time_horzsi = 1:num_time_horzs
        
        %pull out time segments of variable of interest in discount rate
            time_start_ind = 1;
            time_end_ind = find(time_desc == time_horzs(time_horzsi));
        
            per_capita_GDP_time_seg = per_capita_GDP_thousand_per_person(time_start_ind:time_end_ind,:);
            
            R_time_seg = R(discounti,time_start_ind:time_end_ind)';
            
        %discount each time segment (for each mitigation level)
            per_capita_GDP_time_seg_disc = NaN(size(per_capita_GDP_time_seg));
            
            for exp_i = 1:num_scens
                
                per_capita_GDP_time_series_this_exp = squeeze(per_capita_GDP_time_seg(:,exp_i));

                %apply discount
                %per_capita_GDP_time_series_this_exp_disc = (per_capita_GDP_time_series_this_exp.*R_time_seg)./nansum(R_time_seg);
                per_capita_GDP_time_series_this_exp_disc = (per_capita_GDP_time_series_this_exp.*R_time_seg);
                
                %store discounted time series
                per_capita_GDP_time_seg_disc(:,exp_i) = per_capita_GDP_time_series_this_exp_disc;

            end
            
            %integrate over time

            npv_per_capita_GDP_time_seg_disc = nansum(per_capita_GDP_time_seg_disc,1);
            
            %find mitigation level with optimim GDP per capita

            ind_with_max_GDP_per_cap = find(npv_per_capita_GDP_time_seg_disc == max(npv_per_capita_GDP_time_seg_disc));
            
            %put in the optimum mitigation level (labelled with max temp)
        
            optimim_level_of_mitigation_stab_temp(discounti,time_horzsi) = max_T(ind_with_max_GDP_per_cap(1));
        
    end
end

%% plot

two_d_plot = optimim_level_of_mitigation_stab_temp;
x_axis = time_horzs;
y_axis = discount_rates;

FigHandle = figure('Position', [100, 100, 700, 600]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName', 'helvetica')
hold on
grid on
%colormap(redblue)
%colormap(redgreencmap)
%colormap(pink)

                max_color_range = 4.0;
                min_color_range = 1.3;
                num_contours = 20;

                contourf(x_axis,y_axis,two_d_plot,linspace(min_color_range,max_color_range,num_contours), 'LineStyle','none');
                contour(x_axis,y_axis,two_d_plot,[1.5 2 3],'k','LineStyle','-','LineWidth',4);
                caxis([min_color_range max_color_range])
                
                plot([2100 2100],[min(discount_rates) max(discount_rates)],'--k')
                plot([2300 2300],[min(discount_rates) max(discount_rates)],'--k')

                plot([min(time_horzs) max(time_horzs)],[0.0 0.0],'--k')
                plot([min(time_horzs) max(time_horzs)],[0.03 0.03],'--k')
                %plot([min(time_horzs) max(time_horzs)],[0.05 0.05],'--k')
                %plot([min(time_horzs) max(time_horzs)],[0.07 0.07],'--k')
                
                %sanePColor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');
                %pcolor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');

                h = colorbar;

                xlim([min(time_horzs) max(time_horzs)])
                ylim([min(discount_rates) max(discount_rates)])
                
%                 h.Ruler.Scale = 'log';
%                 h.Ruler.MinorTick = 'on';

                ylabel(h,'Optimal level of mitigation effort (global warming, C)');
                
                ylabel('Discount rate')
                xlabel('Time horizon (year)')
                
                %title('optimal mitigation level f(discount rate,time horizon)')

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/opt_level_mitigation_time_horz_discount','-dpdf')