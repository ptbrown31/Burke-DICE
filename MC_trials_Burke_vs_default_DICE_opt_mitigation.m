close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

%% make set MC paramaters

    T = 65;
    
    range_values_2nd_half = linspace(1,1.5,6);
    range_values_1st_half = 1./(range_values_2nd_half(2:end));

    range_values = horzcat(fliplr(range_values_1st_half),range_values_2nd_half);
    
    num_mc_trials = 2000;
    num_varying_variables = 8;
    
   mc_var_names = {'CO2 intensity of economy, growth rate',...
                   'CO2 intensity of economy, \Delta growth rate',...
                   'backstop cost decline rate',...
                   'mitigation cost curve exponent',...
                   'asymptotic population',...
                   'TFP, growth rate',...
                   'TFP, \Delta growth rate',...
                   '2XCO2 climate sensitivity'};
    
    index_pull = randi(length(range_values),num_varying_variables,num_mc_trials);  %note that index pul is same for both damage functions
    
    %generate the indicies pulled from each varying variable

overall_plot_series_all_all = NaN(17,T,61,num_mc_trials,2);


lin_count = 1;

for damage_i = 0:1
    for mc_trial_i = 1:num_mc_trials
        
        disp('fraction done')
        lin_count./(2*num_mc_trials)

    global abatement_cost_on
    abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 

    global Burke_damage_on;
    Burke_damage_on = damage_i;  %1 means Burke damage through TFP and capital depreciation pathways, 
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

%revove the 1st 4 that are too fast
aopts = aopts(:,5:end);

num_scens = size(aopts,2);

    %% run main program

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
    
    middle_value = 6;

                               damage_dep_rate_coeff_on_T = damage_dep_rate_coeffs_on_T(middle_value);
                               damage_TFP_growth_coeff_on_T = damage_TFP_growth_coeffs_on_T(middle_value);
                               elasmu = elasmu_vals(middle_value);
                               Sigg0 = Sigg0_vals(index_pull(1,mc_trial_i));
                               Siggd = Siggd_vals(index_pull(2,mc_trial_i));
                               pd = pd_vals(index_pull(3,mc_trial_i));      
                               theta2 = theta2_vals(index_pull(4,mc_trial_i));
                               LA = LA_vals(index_pull(5,mc_trial_i));
                               Ag0 = Ag0_vals(index_pull(6,mc_trial_i));
                               Agd = Agd_vals(index_pull(7,mc_trial_i));
                               prstp = prstp_vals(middle_value);
                               deltaT = deltaT_vals(index_pull(8,mc_trial_i));
                               

                               overall_plot_series_all = [];

                               for aopti = 1:num_scens
                                 
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

    % Net output: 9 ($ trillion 2005 USD)
    % Population: 7 (million people)
    % Output/capita = $ million/person

    per_capita_GDP_million_per_person = overall_plot_series_all(9,:,:)./overall_plot_series_all(7,:,:);

    %convert to thousand / person

    per_capita_GDP_thousand_per_person = (1000).*per_capita_GDP_million_per_person;

    overall_plot_series_all = cat(1,overall_plot_series_all,per_capita_GDP_thousand_per_person);
    
    overall_plot_series_all_all(:,:,:,mc_trial_i,damage_i+1) = overall_plot_series_all;

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
              'fract output w.r.t. no abatement',...
              '$ thousand 2005 USD / person'};
          
          lin_count = lin_count + 1;
          
    end   
end

%% save

value_of_varying_variables = cat(1, Sigg0_vals,...
                                    Siggd_vals,...
                                    pd_vals,...
                                    theta2_vals,...
                                    LA_vals,...
                                    Ag0_vals,...
                                    Agd_vals,...
                                    deltaT_vals);

save('/home/pbrown/Documents/MATLAB/DICE/test_Burke/MC_trials_Burke_vs_default_DICE_opt_mitigation.mat',...
      'overall_plot_series_all_all',...
      'titles',...
      'ylabels',...
      'index_pull',...
      'value_of_varying_variables',...
      'mc_var_names',...
      '-v7.3')

