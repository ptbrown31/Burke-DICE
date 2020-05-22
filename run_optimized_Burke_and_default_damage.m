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

global optimize_on;
optimize_on = 1;  %1 means optimizing, 0 means its reading in a control rate

%     if optimize_on == 0;
%         
         global aopt
%         aopt = ones(T,1);
%         
%     end

global optimize_for_total_utility;
optimize_for_total_utility = 1; %1 means optimize for total utility, 0 means optimize for per-capita utility

global utility_function;
utility_function = 1; %1 means utility function in code and 0 means utility function in text,

%% make a bunch of aopts to feed the model

time = 1:T;
tstep = 5; %  (5 years per period)

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
Agd_vals = range_values*0.005;           % 10

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

                                %run DICE
                                addpath(genpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth'))

                                DICE;

                                NPV_utility = fval;

                                %Calculate the Carbon Price
                                opt_price = pbacktime.*(aopt'./pai(1:length(aopt))).^(theta2-1); % Carbon price, $ per ton of CO2

                                overall_plot_series_all = cat(1,S(16,:),...
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
      
% Net output: 9 ($ trillion 2005 USD)
% Population: 7 (million people)
% Output/capita = $ million/person
   
per_capita_GDP_million_per_person = overall_plot_series_all(9,:,:)./overall_plot_series_all(7,:,:);

%convert to thousand / person

per_capita_GDP_thousand_per_person = (1000).*per_capita_GDP_million_per_person;

overall_plot_series_all = cat(1,overall_plot_series_all,per_capita_GDP_thousand_per_person);

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
          'abatement cost fraction',...
          'per-capita GWP'};
      
ylabels = {'utils',...
          'utils / person',...
          '$ trillion 2010 USD',...
          '$ thousand 2010 USD / person',...
          '$ trillion 2010 USD',...
          'MTC/$1000',...
          'million people',...
          'fraction',...
          '$ trillion 2010 USD',...
          'Gt CO_2 / year',...
          'Gt CO_2',...
          '\Delta T above PI',...
          'frac output w.r.t. no warming',...
          'fraction',...
          '$1000 / Ton CO_2',...        
          'fract output w.r.t. no abatement',...
          '$ thousand 2010 USD / person'};
      


%% save stuff for other programs

    % find the aopt index that is optimum
    
    aopt_time_series = overall_plot_series_all(14,:);

code_string = ['_a',num2str(abatement_cost_on),...
                     '_B',num2str(Burke_damage_on),...
                     '_o',num2str(optimize_on),...
                     '_u',num2str(optimize_for_total_utility),...
                     '_uf',num2str(utility_function)];
                 
% strcat('overall_plot_series_',code_string) = overall_plot_series;
% strcat('NPV_utility_',code_string) = NPV_utility;

      save(['/home/pbrown/Documents/MATLAB/DICE/test_Burke/scale_with_temperature_Burke_damage_aopt',code_string,'.mat'],...
          'aopt_time_series',...
          'overall_plot_series_all',...
          'titles',...
          'ylabels',...
          'code_string')

