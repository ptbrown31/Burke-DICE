close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

%% decide on type of experiment 

global abatement_cost_on
abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 

global Burke_damage_on;
Burke_damage_on = 1;  %1 means Burke damage through TFP and capital depreciation pathways, 
                      %0 means default DICE, 
                      %2 means no damage
                      %3 means damage with a threshold for 2C so it won't go
                      %over that value

global optimize_on;
optimize_on = 1;  %1 means optimizing, 
                  %0 means its reading in a control rate for 1.5C
                  %2 means its a control rate for 2.0C

    if optimize_on == 0;
        
        load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard_a0_B1_o1_u1.mat',...
        'overall_plot_series')
    
        global aopt

        aopt = overall_plot_series(14,:)';

    end
    
    if optimize_on == 2;
        
        load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard_a1_B1_o1_u1.mat',...
        'overall_plot_series')
    
        global aopt

        aopt = overall_plot_series(14,:)';
        
        aopt = 1.4*aopt;
        
        aopt(aopt>1) = 1;

    end
    
global optimize_for_total_utility;
optimize_for_total_utility = 1; %1 means optimize for total utility, 0 means optimize for per-capita utility

global utility_function;
utility_function = 1; %1 means utility function in code and 0 means utility function in text, 2 means utility linear with per-capita consumption

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
prstp_vals = range_values*0.015;         % 11

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
    
               
%% plot 

time = 2010:5:2330;
      
row_num = 4;
col_num = 4;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',10);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
    end_year_plot = 2250;

for plot_i = 1:length(titles)
    
    subplot(row_num,col_num,plot_i)
    hold on
    grid on
    
        %plot all
        
        %plot highlights
        plot(time,squeeze(overall_plot_series(plot_i,:)),'-k','LineWidth',1.5)
              
        xlim([2010 end_year_plot])

        xlabel('Year')
        ylabel(ylabels{plot_i})
        title(titles{plot_i})
        
end  

%% save these results

code_string = ['_a',num2str(abatement_cost_on),...
                     '_B',num2str(Burke_damage_on),...
                     '_o',num2str(optimize_on),...
                     '_u',num2str(optimize_for_total_utility),...
                     '_uf',num2str(utility_function)];
                 
% strcat('overall_plot_series_',code_string) = overall_plot_series;
% strcat('NPV_utility_',code_string) = NPV_utility;

      save(['/home/pbrown/Documents/MATLAB/DICE/test_Burke/Burke_damage_series_standard',code_string,'.mat'],...
          'overall_plot_series',...
          'NPV_utility',...
          'code_string')