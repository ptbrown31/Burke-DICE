close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/')

%% read in Burke data
row_offset = 0;
col_offset = 0;

rcps = {'RCP26','RCP45','RCP60','RCP85'};
time_horzs = {'2049','2099'};

filename_base = '/home/pbrown/Documents/MATLAB/DICE/test_Burke/BDD18_Fig_4a_data/BDD18_Fig4a_GDP_per_cap_vs_GMSAT_SSP1_';

GDP_loss_Burke_Fig_4a_all_data = NaN(39,2,4,2); %1) number of points (models), 2) GMSAT vs. per capita GDP loss, 3)RCPs, 4) time integration

for time_horizs_i = 1:length(time_horzs)
    for rcp_i = 1:length(rcps)
        
        filename = strcat(filename_base,rcps{rcp_i},'_',time_horzs{time_horizs_i},'.csv');
        
        M = csvread(filename,row_offset,col_offset);
        
        num_rows = size(M,1);
        
        GDP_loss_Burke_Fig_4a_all_data(1:num_rows,:,rcp_i,time_horizs_i) = M;
        
    end
end

zero_inds = find(GDP_loss_Burke_Fig_4a_all_data == 0);

GDP_loss_Burke_Fig_4a_all_data(zero_inds) = NaN;

%% read in SSP1 data 

init_time = horzcat(2005,2010:10:2100);

%IAM: IMAGE

ssp_sub_scen = {'SSP1-26',...
                'SSP1-34',...
                'SSP1-45',...
                'SSP1-Baseline'};
            
var_names = {'GDP',...
             'population',...
             'temperature',...
             'carbon price',...
             'CO_2 emissions'};
         
var_units = {'billion US$2005/yr',...
             'million',...
             'K',...
             'US$2005/t CO2',...
             'Mt CO_2/yr'};
         
row_offset = 0;
col_offset = 0;

filename = '/home/pbrown/Documents/MATLAB/DICE/test_Burke/BDD18_Fig_4a_data/BDD18_Fig4a_SSP1_through_5_data.xlsx';

SSP1_data_raw = NaN(length(ssp_sub_scen),length(init_time),length(var_names)); %1) SSP1 sub scenarios (RCPs), 2) time, 3) variable

%read
    SSP1_data_raw(:,:,1) = xlsread(filename,'SSP1','E5:O8');
    SSP1_data_raw(:,:,2) = xlsread(filename,'SSP1','E11:O14');
    SSP1_data_raw(:,:,3) = xlsread(filename,'SSP1','E17:O20');
    SSP1_data_raw(:,:,4) = xlsread(filename,'SSP1','E23:O26');
    SSP1_data_raw(:,:,5) = xlsread(filename,'SSP1','E29:O32');

    SSP1_data_raw(SSP1_data_raw == -99999) = NaN;
    
%interpolate in time

time_desc = 2005:5:2100;

SSP1_data = NaN(length(ssp_sub_scen),length(time_desc),length(var_names));

for ssp_sub_scen_i = 1:length(ssp_sub_scen)
    for var_names_i = 1:length(var_names)
        
        SSP1_data(ssp_sub_scen_i,:,var_names_i) = interp1(init_time,...
                                                  squeeze(SSP1_data_raw(ssp_sub_scen_i,:,var_names_i)),...
                                                  time_desc);
    end
end

%make a gdp per capita 

gdp_per_capita = SSP1_data(:,:,1)./SSP1_data(:,:,2);

SSP1_data = cat(3,SSP1_data,gdp_per_capita);

var_names = horzcat(var_names,'GDP per capita');

%convert units to same as DICE

% var_units = {'billion US$2005/yr',...
%              'million',...
%              'K',...
%              'US$2005/t CO2',...
%              'Mt CO_2/yr',...
%              'thousand US$2005/yr'};

conversion_farctors = [1/1000 1 1 1 1/1000 1];

var_units = {'trillion US$2005/yr',...
             'million',...
             'K',...
             'US$2005/t CO2',...
             'Gt CO_2/yr',...
             'thousand US$2005/person/yr'};
         
for ssp_sub_scen_i = 1:length(ssp_sub_scen)
    for time_i = 1:length(time_desc)
        
        SSP1_data(ssp_sub_scen_i,time_i,:) = conversion_farctors'.*squeeze(SSP1_data(ssp_sub_scen_i,time_i,:));
        
    end
end         

%% calc production

% OUTPUT AND CAPITAL ACCUMULATION
%   capital share
gama = 0.3;

% PRODUCTIVITY

    T = 65;
    time = 1:T;
    tstep = 5; %  (5 years per period)

    %   Initial level of total factor productivity
    A0 = 3.8; % DICE-2013R; 27.22 for DICE-2007
    %   Initial growth rate for TFP per 5 years
    Ag0 = 0.079; % DICE-2013R
    %   Decline rate of TFP per 5 years
    Agd = 0.006; % DICE-2013R
    %   Rate of growth of productivity (percent per decade)
    Ag = zeros(1,T);
    %   Total factor productivity
    A = zeros(1,T);
    A(1) = A0;
    for i = 1:(T-1)
        Ag(i)=Ag0*exp(-Agd*tstep*(time(i)-1));
        A(i+1) = A(i)/(1-Ag(i));
    end
    
    A = A(1:length(time_desc));
    Ag = Ag(1:length(time_desc));

% have TFP, population, gross output, need to back out kapital

ssp_sub_scen_i = 1;

%http://www.wolframalpha.com/input/?i=solve+for+K,+Y%3DA*(K%5Eg)*(L%5E(1-g))

Y = squeeze(SSP1_data(ssp_sub_scen_i,:,1));
L = squeeze(SSP1_data(ssp_sub_scen_i,:,2));

% this is what the production function looks like: S2(7) = A_damaged * (S2(1) ^ gama) * (L(t1 + 1)/1000)^(1 - gama);

K = ((Y.*((L./1000).^(gama-1)))./(A)).^(1./gama);

% make sure this works by calculating output again

    % Y2 = A .* (K .^ gama) .* (L./1000).^(1 - gama);
    % 
    % Y-Y2
    % it does
    
% calculate what this implies in terms of investiment (savings rate)
    tstep = 5;
    cap_dep_rate = 0.1;

% based on this equation
% S2(1) = tstep * optlrsav * Q + (1 - deltak_damaged) ^ tstep * S1(1);

% K_next = tstep * saverate_now * Y_now + (1 - cap_dep_rate) ^ step * K_now
    % solve for saverate_now ->
    % saverate_now = (K_next - ((1 - cap_dep_rate) ^ step * K_now))./(tstep * Y_now)

    %assume that capital depreciation rate is 0.1
    save_rate = NaN(1,length(time_desc));

    for time_i = 1:length(time_desc)-1
        
        K_now = K(time_i);
        K_next = K(time_i+1);
        Y_now = Y(time_i);
        
        saverate_now = (K_next - ((1 - cap_dep_rate) ^ tstep * K_now))./(tstep * Y_now);
        
        save_rate(time_i) = saverate_now;

    end

time_desc_2050_i = find(time_desc == 2050);
time_desc_2100_i = find(time_desc == 2100);

Y_per_L_2050_no_damage = Y(time_desc_2050_i)./L(time_desc_2050_i);
Y_per_L_2100_no_damage = Y(time_desc_2100_i)./L(time_desc_2100_i);

%% plot these

line_plots = vertcat(Y,A,K,L);

titles = {'gross output',...
          'TFP',...
          'Kapital',...
          'Population'};
      
ylabels = {'trillion US$2005/yr',...
            '-',...
            'trillion US$2005',...
            'million'};

FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',11);
    set(0,'defaultAxesFontName', 'helvetica')
    
for plot_i = 1:length(titles)
    
    subplot(2,2,plot_i)
    hold on
    grid on
    
    plot(time_desc,line_plots(plot_i,:),'-k','LineWidth',2)
    
    title(titles{plot_i})
    
    ylabel(ylabels{plot_i})
    xlabel('year')

end

%% tune it

range_values_2nd_half = linspace(1,1.5,6);
range_values_1st_half = 1./(range_values_2nd_half(2:end));

range_values = horzcat(fliplr(range_values_1st_half),range_values_2nd_half);

%choose a central coeff

% damage_dep_rate_coeffs_on_T = range_values*0.063;
% damage_TFP_growth_coeffs_on_T = range_values*0.0025;

    damage_dep_rate_coeffs_on_T = range_values*(1)*0.063;
    damage_TFP_growth_coeffs_on_T = range_values*(1)*0.0025;

    middle_value_i = 6; 

    damage_dep_rate_coeff_on_T = damage_dep_rate_coeffs_on_T(middle_value_i);
    damage_TFP_growth_coeff_on_T = damage_TFP_growth_coeffs_on_T(middle_value_i);

%I have 4 global temperature trajectories
%this are connected to TFP and K depreciation

Y_per_Ls_vs_no_damage = NaN(length(ssp_sub_scen),length(time_desc));

%need to run model forward in time such that it produces correct losses

deltak = 0.1;
optlrsav = 0.258278145695364;

A_with_damage = NaN(length(ssp_sub_scen),length(time_desc));
K_with_damage = NaN(length(ssp_sub_scen),length(time_desc));
Y_with_damage = NaN(length(ssp_sub_scen),length(time_desc));
Y_per_L_with_damage = NaN(length(ssp_sub_scen),length(time_desc));

Y_without_damage = NaN(length(ssp_sub_scen),length(time_desc));
Y_per_L_without_damage = NaN(length(ssp_sub_scen),length(time_desc));

A_with_damage(:,1) = A(1);
K_with_damage(:,1) = K(1);
Y_with_damage(:,1) = Y(1);
Y_per_L_with_damage(:,1) = Y(1)./squeeze(SSP1_data_raw(:,1,2));

Y_per_L_without_damage(:,1) = Y(1)./squeeze(SSP1_data_raw(:,1,2));
Y_without_damage(:,1) = Y(1)./squeeze(SSP1_data_raw(:,1,2));

% var_names = {'GDP',...
%              'population',...
%              'temperature',...
%              'carbon price',...
%              'CO_2 emissions'};

% SSP1_data_raw = NaN(length(ssp_sub_scen),length(init_time),length(var_names)); %1) SSP1 sub scenarios (RCPs), 2) time, 3) variable

for ssp_sub_scen_i = 1:length(ssp_sub_scen)
    
    %define intermediate things you dont need to save
    %Ag_damaged = NaN(1,length(time_desc));
    %A_damaged = NaN(1,length(time_desc));    
        
    %deltak_damaged = NaN(1,length(time_desc)-1);
    
    for time_i = 2:length(time_desc)
        
        %what is the global warming?
        
        global_temp_now = SSP1_data(ssp_sub_scen_i,time_i,3);
        
        %calculate the damaged TFP for this temperature

%                Ag_damaged = Ag(t1-1) - damage_TFP_growth_coeff_on_T*S1(2);
%                A_damaged = S1(23)/(1-Ag_damaged);  %S1(23) is the previous timesteps damaged TFP

                 Ag_damaged_now = Ag(time_i-1) - damage_TFP_growth_coeff_on_T*global_temp_now;
                 
                 A_with_damage(ssp_sub_scen_i,time_i) = A_with_damage(ssp_sub_scen_i,time_i-1)./(1-Ag_damaged_now);
%                 A_with_damage(ssp_sub_scen_i,time_i) = A(time_i);
        
        %calculate the damaged deltak for this temperature
        
%               deltak_damaged = damage_dep_rate_coeff_on_T*S1(2);
%               if deltak_damaged < deltak; deltak_damaged = deltak; end

                deltak_damaged_now = damage_dep_rate_coeff_on_T*global_temp_now;
                if deltak_damaged_now < deltak; deltak_damaged_now = deltak; end

        %calculate the damaged capital stock (dont use optlrsave but
        %emperically derived time-varying saveings rate
        
%               S2(1) = tstep * optlrsav * Q + (1 - deltak_damaged) ^ tstep * S1(1);
 
%                 K_with_damage(ssp_sub_scen_i,time_i) = tstep * optlrsav * Y_with_damage(ssp_sub_scen_i,time_i-1) + ...
%                                                        (1 - deltak_damaged_now) ^ tstep * K_with_damage(ssp_sub_scen_i,time_i-1);

                K_with_damage(ssp_sub_scen_i,time_i) = tstep * save_rate(time_i-1) * Y_with_damage(ssp_sub_scen_i,time_i-1) + ...
                                                       (1 - deltak_damaged_now) ^ tstep * K_with_damage(ssp_sub_scen_i,time_i-1);

        %calculate the damaged gross output
                
%               S2(7) = A_damaged * (S2(1) ^ gama) * (L(t1 + 1)/1000)^(1 - gama);
            
                Y_with_damage(ssp_sub_scen_i,time_i) = A_with_damage(ssp_sub_scen_i,time_i) * ...
                                                      (K_with_damage(ssp_sub_scen_i,time_i) ^ gama) * ...
                                                      (SSP1_data(ssp_sub_scen_i,time_i,2)/1000)^(1 - gama);
                                                  
        %calculate undamaged gross output
        
                Y_without_damage(ssp_sub_scen_i,time_i) = A(time_i) * ...
                                                          (K(time_i) ^ gama) * ...
                                                          (SSP1_data(ssp_sub_scen_i,time_i,2)/1000)^(1 - gama);
            
        %calculate the damaged gross output per capita
        
                Y_per_L_with_damage(ssp_sub_scen_i,time_i) = Y_with_damage(ssp_sub_scen_i,time_i)./SSP1_data(ssp_sub_scen_i,time_i,2);
                
                Y_per_L_without_damage(ssp_sub_scen_i,time_i) = Y_without_damage(ssp_sub_scen_i,time_i)./SSP1_data(ssp_sub_scen_i,time_i,2);
        
                Y_per_Ls_vs_no_damage(ssp_sub_scen_i,time_i) = -100*(1-(Y_per_L_with_damage(ssp_sub_scen_i,time_i)./Y_per_L_without_damage(ssp_sub_scen_i,time_i)));
        
    end
end

Y_per_Ls_vs_no_damage_points = Y_per_Ls_vs_no_damage(:,[time_desc_2050_i time_desc_2100_i]);
Temp_points = SSP1_data(:,[time_desc_2050_i time_desc_2100_i],3);

%% plot the components of where the damage is coming from 

FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',11);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
    subplot(2,2,1)
    hold on
    grid on
    
        plot(time_desc,A,'-k','LineWidth',2)

        for rcp_i = 1:length(rcps)
            plot(time_desc,A_with_damage(rcp_i,:),'-r','LineWidth',1)
        end

        xlabel('Year')
        ylabel('X')
        title('A: TFP')
        
    subplot(2,2,2)
    hold on
    grid on
    
        plot(time_desc,K,'-k','LineWidth',2)

        for rcp_i = 1:length(rcps)
            plot(time_desc,K_with_damage(rcp_i,:),'-r','LineWidth',1)
        end

        xlabel('Year')
        ylabel('X')
        title('K: Capital')
        
    subplot(2,2,3)
    hold on
    grid on
    
        plot(time_desc,Y,'-k','LineWidth',2)

        for rcp_i = 1:length(rcps)
            plot(time_desc,Y_with_damage(rcp_i,:),'-r','LineWidth',1)
        end

        xlabel('Year')
        ylabel('X')
        title('Y: Output')
        
    subplot(2,2,4)
    hold on
    grid on
    
        plot(time_desc,(Y./L),'-k','LineWidth',2)

        for rcp_i = 1:length(rcps)
            plot(time_desc,Y_per_L_with_damage(rcp_i,:),'-r','LineWidth',1)
        end

        xlabel('Year')
        ylabel('X')
        title('Y/L: Per-Capita Output')         
%% plot

colors = {'g','r','b','m'};
syms = {'o','o'};

FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',11);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
for time_horizs_i = 1:length(time_horzs)
    for rcp_i = 1:length(rcps)
        if time_horizs_i == 1;
            scatter(GDP_loss_Burke_Fig_4a_all_data(:,1,rcp_i,time_horizs_i),...
                    GDP_loss_Burke_Fig_4a_all_data(:,2,rcp_i,time_horizs_i),100,colors{rcp_i},syms{time_horizs_i});
        end
        if time_horizs_i == 2;
            scatter(GDP_loss_Burke_Fig_4a_all_data(:,1,rcp_i,time_horizs_i),...
                    GDP_loss_Burke_Fig_4a_all_data(:,2,rcp_i,time_horizs_i),100,colors{rcp_i},syms{time_horizs_i},'filled');
        end
        
        scatter(Temp_points(rcp_i,time_horizs_i),Y_per_Ls_vs_no_damage_points(rcp_i,time_horizs_i),200,'k','s','filled')
        
    end
end

ylabel('Change in global GDP (%)')
xlabel('Global warming (C)')

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/tune_Burke_damage_paramaters','-dpdf')
















%% decide on type of experiment 

% T = 65;
% 
% global abatement_cost_on
% abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 
% 
% global Burke_damage_on;
% Burke_damage_on = 1;  %1 means Burke damage through TFP and capital depreciation pathways, 
%                       %0 means default DICE, 
%                       %2 means no damage
%                       %3 means damage with a threshold for 2C so it won't go
%                       %over that value
%                       
% global optimize_on;
% optimize_on = 0;  %1 means optimizing, 0 means its reading in a control rate
% 
% %     if optimize_on == 0;
% %         
%          global aopt
% %         aopt = ones(T,1);
% %         
% %     end
% 
% global optimize_for_total_utility;
% optimize_for_total_utility = 1; %1 means optimize for total utility, 0 means optimize for per-capita utility
% 
% global utility_function;
% utility_function = 1; %1 means utility function in code and 0 means utility function in text, 2 means utility linear with per-capita consumption
% 
% %% make a bunch of aopts to feed the model
% 
% time = 1:T;
% tstep = 5; %  (5 years per period)
% 
% %make as many control rates as there are timesteps: 1 reaching 1 in each
% %timestep
% 
% aopts = ones(T,T);
% 
% for aopti = 1:T
%    
%     aopts(1:aopti,aopti) = linspace(0,1,aopti);
%     
% end
% 
% aopts(:,end) = zeros(T,1); %make the last one no climate policy
% 
% 
% %% run main program
% 
% range_values_2nd_half = linspace(1,1.5,6);
% range_values_1st_half = 1./(range_values_2nd_half(2:end));
% 
% range_values = horzcat(fliplr(range_values_1st_half),range_values_2nd_half);
% 
% %accordian the variables
% 
% %damage function
% damage_dep_rate_coeffs_on_T = range_values*0.063;
% damage_TFP_growth_coeffs_on_T = range_values*0.0025;
% 
% %Elasticity of marginal utility 
% %of consumption
% elasmu_vals = range_values*1.45;         % 3
% elasmu_vals(elasmu_vals <= 1) = 1.001;   
% 
% % carbon intensity of output
% Sigg0_vals = range_values*-0.0152;         % 4
% Siggd_vals = range_values*-0.001;        % 5
% 
% %initial cost decline backstop 
% %cost per period
% pd_vals = range_values*0.025;            % 6
% 
% %abatement cost curve exponent
% theta2_vals = range_values*2.6;          % 7
% 
% %asymptotic population
% LA_vals = range_values*11500;            % 8
% 
% %total factor productivity
% Ag0_vals = range_values*0.079;           % 9
% Agd_vals = range_values*0.006;           % 10
% 
% %pure rate of time preference
% prstp_vals = range_values*0.03;         % 11
% 
% %climate sensitivity
% deltaT_vals = range_values*3.1;          % 12
% 
% num_vars = 12;
% 
% var_names = {'damage cap-dep-rate coeff on T',...
%              'damage TFP coeff on T',...
%              'elasticity of utility',...
%              'C intensity, growth rate',...
%              'C intensity, change in growth rate',...
%              'backstop cost decline',...
%              'abatement cost curve exponent',...
%              'asymptotic population',...
%              'TFP, growth rate',...
%              'TFP, change in growth rate',...
%              'pure rate of time preference',...
%              'climate sensitivity'};
% 
% global damage_dep_rate_coeff_on_T
% global damage_TFP_growth_coeff_on_T
% global elasmu
% global Sigg0
% global Siggd
% global pd
% global theta2
% global LA
% % global Ag0
% % global Agd
% global prstp
% global deltaT
% 
% middle_value_i = 6;
% 
%         %set everything to middle
%         
%                            damage_dep_rate_coeff_on_T = damage_dep_rate_coeffs_on_T(middle_value_i);
%                            damage_TFP_growth_coeff_on_T = damage_TFP_growth_coeffs_on_T(middle_value_i);
%                            elasmu = elasmu_vals(middle_value_i);
%                            Sigg0 = Sigg0_vals(middle_value_i);
%                            Siggd = Siggd_vals(middle_value_i);
%                            pd = pd_vals(middle_value_i);      
%                            theta2 = theta2_vals(middle_value_i);
%                            LA = LA_vals(middle_value_i);    
% %                            Ag0 = Ag0_vals(middle_value_i);
% %                            Agd = Agd_vals(middle_value_i);    
%                            prstp = prstp_vals(middle_value_i);
%                            deltaT = deltaT_vals(middle_value_i);
%                            
%                            overall_plot_series_all = [];
%                            
%                            for aopti = 1:T
%                                aopti
%                                
%                                aopt = aopts(:,aopti);
% 
%                                 %run DICE
%                                 addpath(genpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth_tune_parameters'))
% 
%                                 DICE;
% 
%                                 NPV_utility = fval;
% 
%                                 %Calculate the Carbon Price
%                                 opt_price = pbacktime.*(aopt'./pai(1:length(aopt))).^(theta2-1); % Carbon price, $ per ton of CO2
% 
%                                 overall_plot_series = cat(1,S(16,:),...
%                                                     S(17,:),...
%                                                     S(7,:),...
%                                                     S(19,:),...
%                                                     S(9,:),...
%                                                     Sig,...
%                                                     L,...
%                                                     R,...
%                                                     S(20,:),...
%                                                     S(12,:),...
%                                                     S(4,:),...
%                                                     S(2,:),...
%                                                     S(14,:),...
%                                                     aopt',...
%                                                     opt_price,...
%                                                     S(10,:));
%                                                 
%                                                 overall_plot_series_all = cat(3,overall_plot_series_all,overall_plot_series);
%                                                 
%                            end
%     % S
%     % 1 - regular capital X
%     % 2 - GMST X
%     % 3 - lower ocean temperature X
%     % 4 - atmospheric CO2 X
%     % 5 - upper ocean CO2 X
%     % 6 - lower ocean CO2 X
%     % 7 - gross output X
%     % 8 - zero-carbon emissions captial (leave NaN)
%     % 9 - total capital (same as 1)
%     % 10 - abatement cost
%     % 11 - epsilon (leave NaN)
%     % 12 - Industrial CO2 emissions
%     % 13 - CO2 Radiative Forcing
%     % 14 - Climate Damage
%     % 15 - 
%     % 16 - total utility
%     % 17 - per capita utility
%     % 18 - consumption
%     % 19 - consumption per capita
%     % 20 - net output
%     % 21 - deltaK
%     % 22 - TFP growth
%     % 23 - TFP
%     
% titles = {'total utility',...
%           'per-capita utility',...
%           'gross output',...
%           'per-capita consumption',...
%           'total capital',...
%           'CO_2 intensity of output',...
%           'population',...
%           'discount factor',...
%           'net output (Y_G)',...
%           'industrial CO_2 emissions',...
%           'ATM CO_2',...
%           'global temperature',...
%           'damage from warming fraction',...
%           'control rate',...
%           'carbon tax',...
%           'abatement cost fraction'};
%       
% ylabels = {'utils',...
%           'utils / person',...
%           '$ trillion 2005 USD',...
%           '$ thousand 2005 USD / person',...
%           '$ trillion 2005 USD',...
%           'MTC/$1000',...
%           'million people',...
%           'fraction',...
%           '$ trillion 2005 USD',...
%           'Gt CO_2 / year',...
%           'Gt CO_2',...
%           '\Delta T above PI',...
%           'frac output w.r.t. no warming',...
%           'fraction',...
%           '$1000 / Ton CO_2',...        
%           'fract output w.r.t. no abatement'};
%       
% %% make discount rates to use after the fact (discount rate isnt used in the above anyway
% 
% % discount_rates = [0 0.01 0.02 0.03 0.05];
% % 
% % discount_rates_styles_1 = {'-k','-r','-b','-g','-c'};
% % discount_rates_styles_2 = {'--k','--r','--b','--g','--c'};
% 
% discount_rates = [0 0.03 0.05];
% legend_text = {'0%','3%','5%'};
% 
% legend_text_full = {'Through 2100, 0% discount',...
%                     'Through 2100, 3% discount',...
%                     'Through 2100, 5% discount',...
%                     'Through 2150, 0% discount',...
%                     'Through 2150, 3% discount',...
%                     'Through 2150, 5% discount'};
% 
% discount_rates_styles_1 = {'-k','-r','-b',};
% discount_rates_styles_2 = {'--k','--r','--b'};
% 
% discount_rates_colors = {'k','r','b',};
% discount_rates_symbols_1 = {'s','s','s'};
% discount_rates_symbols_2 = {'d','d','d'};
% 
% %   Social rate of time preference (percent per year)
% 
% ro = NaN(length(discount_rates),T);
% 
% for discount_rate_i = 1:length(discount_rates)
% 
%     ro(discount_rate_i,:) = discount_rates(discount_rate_i).*100.*ones(1,T);
% end
% %   Social time preference factor
% R = zeros(length(discount_rates),T);
% R(:,1) = 1;
% for discount_rate_i = 1:length(discount_rates)
%     for i = 2:T
%         R(discount_rate_i,i) = R(discount_rate_i,i-1)/(1+ro(discount_rate_i,i-1)./100).^tstep;
%     end
% end
% 
% %lambda = 1/(1+ro(:,1)/100)^tstep;
%       
% %% get a max temp for each aopt
% 
% time_desc = 2010:5:2330;
% 
% max_T = NaN(T,1);
% 
% %diff everything w.r.t caase [b]    
%     for exp_i = 1:T
% 
%         max_T(exp_i) = max(squeeze(overall_plot_series_all(12,:,exp_i)));
% 
%     end   
% 
%     FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',10);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
%     grid on    
%     
%     plot(time_desc,max_T,'-k','LineWidth',2)
%     scatter(time_desc,max_T,25,[0 0 0],'filled')
%     
%     ylabel('maximum temperature over whole period')    
%     xlabel('year of net zero emission')
% 
% 
% %% plot time series for each temp stabilization
% 
% row_num = 4;
% col_num = 4;
% 
% end_year_plot = 2150;
% 
%     FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',6);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
%     grid on
%     
% for var_i = 1:length(titles)
%     
%     subplot(row_num,col_num,var_i)
%     hold on
%     grid on
%     
%         for exp_i = 1:T
%     
%             %plot highlights
%             plot(time_desc,squeeze(overall_plot_series_all(var_i,:,exp_i)),'LineWidth',0.5)
%         
%         end
%               
%         xlim([2010 end_year_plot])    
% 
%         %xlim([2015 end_year_plot])
% 
%         xlabel('time')
%         ylabel(ylabels{var_i})
%         title(titles{var_i})
% 
% end
% 
% print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/time_series_all_temp_targets','-dpdf')
% 
% %% plot tight
% 
% sub_var_inds = [10 12 16 13 2 4];
% 
% row_num = 3;
% col_num = 2;
% 
% end_year_plot = 2150;
% 
%     FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
%     set(gcf,'color',[1 1 1]);
%     set(0, 'DefaultAxesFontSize',6);
%     set(0,'defaultAxesFontName', 'helvetica')
%     hold on
%     grid on
%     
% for var_i = 1:length(sub_var_inds)
%     
%     subplot(row_num,col_num,var_i)
%     hold on
%     grid on
%     
%         for exp_i = 1:T
%     
%             if var_i ~= 4; plot(time_desc,squeeze(overall_plot_series_all(sub_var_inds(var_i),:,exp_i)),'LineWidth',0.5); end
%             if var_i == 4; plot(time_desc,-1*squeeze(overall_plot_series_all(sub_var_inds(var_i),:,exp_i)-1),'LineWidth',0.5); end
%         
%         end
%               
%         xlim([2010 end_year_plot])    
% 
%         %xlim([2015 end_year_plot])
% 
%         xlabel('time')
%         ylabel(ylabels{sub_var_inds(var_i)})
%         title(titles{sub_var_inds(var_i)})
% 
% end
% 
% print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/time_series_tight_temp_targets','-dpdf')
% 
% %% anom everything w.r.t no policy quickest 100% control and make 2d plot
% 
% overall_plot_series_all_diff = NaN(size(overall_plot_series_all));
% 
%     for exp_i = 1:T
% 
%         overall_plot_series_all_diff(:,:,exp_i) = (squeeze(overall_plot_series_all(:,:,exp_i)) - squeeze(overall_plot_series_all(:,:,end)));
% 
%     end
%     
% var_to_plot_i = 4; %4 is per capita consumption
% 
% two_d_plot = squeeze(overall_plot_series_all_diff(var_to_plot_i,:,:));
% 
% %% make discounted integrations
% 
% end_year_plot = 2160;
% end_year_calc_1 = 2100;
% end_year_calc_2 = 2150;
% 
% end_year_plot_ind = find(time_desc == end_year_plot);
% end_year_calc_ind_1 = find(time_desc == end_year_calc_1);
% end_year_calc_ind_2 = find(time_desc == end_year_calc_2);
% 
% start_time_ind = 1;
% start_temp_ind = 2;
% 
% two_d_plot_trimmed_1 = two_d_plot(start_time_ind:end_year_calc_ind_1,start_temp_ind:end);
% two_d_plot_trimmed_2 = two_d_plot(start_time_ind:end_year_calc_ind_2,start_temp_ind:end);
% 
% two_d_plot_trimmed_discounted_1 = NaN(size(two_d_plot_trimmed_1,1),size(two_d_plot_trimmed_1,2),length(discount_rates));
% two_d_plot_trimmed_discounted_2 = NaN(size(two_d_plot_trimmed_2,1),size(two_d_plot_trimmed_2,2),length(discount_rates));
% 
% for discount_rate_i = 1:length(discount_rates)
%     for exp_i = 1:T-start_temp_ind
% 
%         two_d_plot_trimmed_discounted_1(:,exp_i,discount_rate_i) = squeeze(two_d_plot_trimmed_1(:,exp_i)).*squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_1,1))')...
%             ./nansum(squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_1,1))')); %weigh it
%         
%         two_d_plot_trimmed_discounted_2(:,exp_i,discount_rate_i) = squeeze(two_d_plot_trimmed_2(:,exp_i)).*squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_2,1))')...
%             ./nansum(squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_2,1))')); %weigh it
%     end
% end
% 
% integrated_net_GDP_1 = NaN(length(max_T(start_temp_ind:end)),length(discount_rates));
% integrated_net_GDP_2 = NaN(length(max_T(start_temp_ind:end)),length(discount_rates));
% 
% for discount_rate_i = 1:length(discount_rates)
%     
%     integrated_net_GDP_1(:,discount_rate_i) = squeeze(nansum(two_d_plot_trimmed_discounted_1(:,:,discount_rate_i),1));
%     integrated_net_GDP_2(:,discount_rate_i) = squeeze(nansum(two_d_plot_trimmed_discounted_2(:,:,discount_rate_i),1));
% 
% end
% 
% int_gdp_min = min(horzcat(min(min(integrated_net_GDP_1)),min(min(integrated_net_GDP_2))));
% int_gdp_max = max(horzcat(max(max(integrated_net_GDP_1)),max(max(integrated_net_GDP_2))));
% 
% %% plot
% 
% FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
% set(gcf,'color',[1 1 1]);
% set(0, 'DefaultAxesFontSize',8);
% set(0,'defaultAxesFontName', 'helvetica')
% %colormap(redblue)
% %colormap(redgreencmap)
% %colormap(pink)
% 
% % rgb = [ ...
% %    158     1    66  
% %    213    62    79   
% %    244   109    67   
% %    253   174    97   
% %    254   224   139   
% %    255   255   191   
% %    230   245   152   
% %    171   221   164
% %    102   194   165   
% %     50   136   189   
% %    94    79   162   ] / 255;
% 
% %b = repmat(linspace(0,1,200),20,1);
% %imshow(b,[],'InitialMagnification','fit')
% %colormap(rgb)
% 
%         subplot(1,2,1)
%         hold on
%         %grid on
%         
%         %plot for legend
%         
%                for discount_rate_i = 1:length(discount_rates)
% 
%                     plot(integrated_net_GDP_1(:,discount_rate_i),...
%                         max_T(start_temp_ind:end),...
%                         discount_rates_styles_1{discount_rate_i},'LineWidth',2)
% 
%                end
%                
%                for discount_rate_i = 1:length(discount_rates)
%                
%                     plot(integrated_net_GDP_2(:,discount_rate_i),...
%                         max_T(start_temp_ind:end),...
%                         discount_rates_styles_2{discount_rate_i},'LineWidth',2)
%                     
%                end
%         
%                %legend(legend_text_full,'location','NorthEast')
%                legend(legend_text_full,'location','NorthWest')
%                
%         %plot again with ID of 'best' temperature 
%                        
%         for discount_rate_i = 1:length(discount_rates)
%             
%             plot(integrated_net_GDP_1(:,discount_rate_i),...
%                 max_T(start_temp_ind:end),...
%                 discount_rates_styles_1{discount_rate_i},'LineWidth',2)
%             
%             plot(integrated_net_GDP_2(:,discount_rate_i),...
%                 max_T(start_temp_ind:end),...
%                 discount_rates_styles_2{discount_rate_i},'LineWidth',2)
%             
%         end
%         
%         for discount_rate_i = 1:length(discount_rates)
%             
%             opt_temp_i_1 = find(squeeze(integrated_net_GDP_1(:,discount_rate_i)) == max(squeeze(integrated_net_GDP_1(:,discount_rate_i))));
%             opt_temp_i_2 = find(squeeze(integrated_net_GDP_2(:,discount_rate_i)) == max(squeeze(integrated_net_GDP_2(:,discount_rate_i))));
%             
%             scatter(integrated_net_GDP_1(opt_temp_i_1,discount_rate_i),...
%                     max_T(opt_temp_i_1),...
%                     100,...
%                     discount_rates_colors{discount_rate_i},...
%                     discount_rates_symbols_1{discount_rate_i},...
%                     'filled')
%                 
%             scatter(integrated_net_GDP_2(opt_temp_i_2,discount_rate_i),...
%                     max_T(opt_temp_i_2),...
%                     100,...
%                     discount_rates_colors{discount_rate_i},...
%                     discount_rates_symbols_2{discount_rate_i},...
%                     'filled')
% 
%         end
%         
%         %plot lines
%         plot([0 0],[min(max_T) max(max_T)],'-k')
%         
%         plot([min([min(min(integrated_net_GDP_1)) min(min(integrated_net_GDP_2))]) max([max(max(integrated_net_GDP_1)) max(max(integrated_net_GDP_2))])],[3 3],'-k')
%         plot([min([min(min(integrated_net_GDP_1)) min(min(integrated_net_GDP_2))]) max([max(max(integrated_net_GDP_1)) max(max(integrated_net_GDP_2))])],[2 2],'-k')
%         plot([min([min(min(integrated_net_GDP_1)) min(min(integrated_net_GDP_2))]) max([max(max(integrated_net_GDP_1)) max(max(integrated_net_GDP_2))])],[1.5 1.5],'-k')
%         
%         xlabel('net per-capita consumption gain w.r.t. no mitigation')
%         ylabel('eventual temperature stabilization (K)')
%         
%         ylim([min(max_T) 0.95*max(max_T)])
%         xlim([1.1*int_gdp_min 1.1*int_gdp_max])
%         
%         subplot(1,2,2)
%         hold on
%         grid on
% 
%                 plot_map_vector = reshape(squeeze(two_d_plot(start_time_ind:end,start_temp_ind:end)),[length(time_desc(start_time_ind:end))*length(max_T(start_temp_ind:end)),1]);
%                 max_color_val = prctile(abs(plot_map_vector),95);
% 
%                 max_color_range = 0.8*max(max(abs(two_d_plot_trimmed_2)));
%                 min_color_range = -1*max_color_range;
%                 num_contours = 25;
% 
%                 contourf(time_desc(start_time_ind:end)',max_T(start_temp_ind:end),two_d_plot(start_time_ind:end,start_temp_ind:end)',linspace(min_color_range,max_color_range,num_contours), 'LineStyle','none');
%                 contour(time_desc(start_time_ind:end)',max_T(start_temp_ind:end),two_d_plot(start_time_ind:end,start_temp_ind:end)',[0 0],'k','LineStyle','-','LineWidth',4);
%                 caxis([min_color_range max_color_range])
%                 
%                 plot([end_year_calc_1 end_year_calc_1],[min(max_T) max(max_T)],'--k')
%                 plot([end_year_calc_2 end_year_calc_2],[min(max_T) max(max_T)],'--k')
% 
%                 %sanePColor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');
% 
%                 %pcolor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');
% 
%                 h = colorbar;
% 
%                 xlim([2020 end_year_plot])
%                 ylim([min(max_T) 0.95*max(max_T)])
%                 
% %                 h.Ruler.Scale = 'log';
% %                 h.Ruler.MinorTick = 'on';
% 
%                 ylabel(h,ylabels{var_to_plot_i});
%                 
%                 ylabel('eventual temperature stabilization (K)')
%                 xlabel('year')
%                 
%                 title(titles{var_to_plot_i})
%                 
