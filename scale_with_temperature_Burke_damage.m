close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

%% decide on type of experiment 

T = 65;

global abatement_cost_on
abatement_cost_on = 1;  %1 means there is an abatement cost and 0 means there is not 

global Burke_damage_on;
Burke_damage_on = 0;  %1 means Burke damage through TFP and capital depreciation pathways, 
                      %0 means default DICE, 
                      %2 means no damage

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
utility_function = 1; %1 means utility function in code and 0 means utility function in text,

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

figure
plot(aopts)


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
      
%% make discount rates to use after the fact (discount rate isnt used in the above anyway

% discount_rates = [0 0.01 0.02 0.03 0.05];
% 
% discount_rates_styles_1 = {'-k','-r','-b','-g','-c'};
% discount_rates_styles_2 = {'--k','--r','--b','--g','--c'};

discount_rates = [0 0.03 0.05];
legend_text = {'0%','3%','5%'};

legend_text_full = {'Through 2100, 0% discount',...
                    'Through 2100, 3% discount',...
                    'Through 2100, 5% discount',...
                    'Through 2335, 0% discount',...
                    'Through 2335, 3% discount',...
                    'Through 2335, 5% discount'};

discount_rates_styles_1 = {'-k','-r','-b',};
discount_rates_styles_2 = {'--k','--r','--b'};

discount_rates_colors = {'k','r','b',};
discount_rates_symbols_1 = {'s','s','s'};
discount_rates_symbols_2 = {'d','d','d'};

%   Social rate of time preference (percent per year)

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

%lambda = 1/(1+ro(:,1)/100)^tstep;
      
%% get a max temp for each aopt

time_desc = 2015:5:2335;

max_T = NaN(num_scens,1);

% %option to find max T over whole time period   
%     for exp_i = 1:num_scens
% 
%         max_T(exp_i) = max(squeeze(overall_plot_series_all(12,:,exp_i)));
% 
%     end

%option to use temperature in 2100
good_time_ind = find(time_desc == 2100);
    for exp_i = 1:num_scens

        max_T(exp_i) = overall_plot_series_all(12,good_time_ind,exp_i);

    end    
    

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',10);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on    
    
    plot(time_desc(end-num_scens+1:end),max_T,'-k','LineWidth',2)
    scatter(time_desc(end-num_scens+1:end),max_T,25,[0 0 0],'filled')
    
    ylabel('maximum temperature over whole period')    
    xlabel('year of net zero emission')


%% plot time series for each temp stabilization

row_num = 4;
col_num = 5;

end_year_plot = 2335;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',6);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
for var_i = 1:length(titles)
    
    subplot(row_num,col_num,var_i)
    hold on
    grid on
    
        for exp_i = 1:num_scens
    
            %plot highlights
            plot(time_desc,squeeze(overall_plot_series_all(var_i,:,exp_i)),'LineWidth',0.5)
        
        end
              
        %xlim([2010 end_year_plot])    

        xlim([2015 end_year_plot])

        xlabel('time')
        ylabel(ylabels{var_i})
        title(titles{var_i})

end

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/time_series_all_temp_targets','-dpdf','-bestfit')

%% plot tight

sub_var_inds = [10 12 16 13 9 17];

row_num = 3;
col_num = 2;

end_year_plot = 2150;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',6);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
for var_i = 1:length(sub_var_inds)
    
    subplot(row_num,col_num,var_i)
    hold on
    grid on
    
        for exp_i = 1:num_scens
    
            if var_i ~= 4; plot(time_desc,squeeze(overall_plot_series_all(sub_var_inds(var_i),:,exp_i)),'LineWidth',0.5); end
            if var_i == 4; plot(time_desc,-1*squeeze(overall_plot_series_all(sub_var_inds(var_i),:,exp_i)-1),'LineWidth',0.5); end
        
        end
              
        %xlim([2010 end_year_plot])    

        xlim([2015 end_year_plot])

        xlabel('time')
        ylabel(ylabels{sub_var_inds(var_i)})
        title(titles{sub_var_inds(var_i)})

end

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/time_series_tight_temp_targets','-dpdf','-bestfit')

%% anom everything w.r.t no policy quickest 100% control and make 2d plot

overall_plot_series_all_diff = NaN(size(overall_plot_series_all));

    for exp_i = 1:num_scens

        overall_plot_series_all_diff(:,:,exp_i) = (squeeze(overall_plot_series_all(:,:,exp_i)) - squeeze(overall_plot_series_all(:,:,end)));

    end
    
var_to_plot_i = 17; %17 is per capita GWP, 9 is net output

two_d_plot = squeeze(overall_plot_series_all_diff(var_to_plot_i,:,:));

%% make discounted integrations

end_year_calc_1 = 2100;
end_year_calc_2 = 2335;
end_year_plot = 2150;

end_year_calc_ind_1 = find(time_desc == end_year_calc_1);
end_year_calc_ind_2 = find(time_desc == end_year_calc_2);
end_year_plot_ind = find(time_desc == end_year_plot);

start_time_ind = 1;
start_temp_ind = 1;

two_d_plot_trimmed_1 = two_d_plot(start_time_ind:end_year_calc_ind_1,start_temp_ind:end);
two_d_plot_trimmed_2 = two_d_plot(start_time_ind:end_year_calc_ind_2,start_temp_ind:end);

two_d_plot_trimmed_discounted_1 = NaN(size(two_d_plot_trimmed_1,1),size(two_d_plot_trimmed_1,2),length(discount_rates));
two_d_plot_trimmed_discounted_2 = NaN(size(two_d_plot_trimmed_2,1),size(two_d_plot_trimmed_2,2),length(discount_rates));

for discount_rate_i = 1:length(discount_rates)
    for exp_i = 1:num_scens-start_temp_ind

%         two_d_plot_trimmed_discounted_1(:,exp_i,discount_rate_i) = squeeze(two_d_plot_trimmed_1(:,exp_i)).*squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_1,1))')...
%             ./nansum(squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_1,1))')); %weigh it
%         
%         two_d_plot_trimmed_discounted_2(:,exp_i,discount_rate_i) = squeeze(two_d_plot_trimmed_2(:,exp_i)).*squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_2,1))')...
%             ./nansum(squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_2,1))')); %weigh it

        two_d_plot_trimmed_discounted_1(:,exp_i,discount_rate_i) = squeeze(two_d_plot_trimmed_1(:,exp_i)).*squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_1,1))');
        
        two_d_plot_trimmed_discounted_2(:,exp_i,discount_rate_i) = squeeze(two_d_plot_trimmed_2(:,exp_i)).*squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted_2,1))');


    end
end

integrated_net_GDP_1 = NaN(length(max_T(start_temp_ind:end)),length(discount_rates));
integrated_net_GDP_2 = NaN(length(max_T(start_temp_ind:end)),length(discount_rates));

for discount_rate_i = 1:length(discount_rates)
    
    integrated_net_GDP_1(:,discount_rate_i) = squeeze(nansum(two_d_plot_trimmed_discounted_1(:,:,discount_rate_i),1));
    integrated_net_GDP_2(:,discount_rate_i) = squeeze(nansum(two_d_plot_trimmed_discounted_2(:,:,discount_rate_i),1));

end

%output the 3% 2100 value for 1.5
[closest_val_15 closest_index_15 ] = min(abs(max_T - 1.5));
[closest_val_2 closest_index_2 ] = min(abs(max_T - 2));
[closest_val_3 closest_index_3 ] = min(abs(max_T - 3));

disp('1.5C, GWP w.r.t. no mitigation, through 2100 3% discount')
integrated_net_GDP_1(closest_index_15,2)

disp('2C, GWP w.r.t. no mitigation, through 2100 3% discount')
integrated_net_GDP_1(closest_index_2,2)

disp('3C, GWP w.r.t. no mitigation, through 2100 3% discount')
integrated_net_GDP_1(closest_index_3,2)

disp('1.5 - 2.0')
integrated_net_GDP_1(closest_index_15,2) - integrated_net_GDP_1(closest_index_2,2)

disp('1.5C, GWP w.r.t. no mitigation, through 2335 3% discount')
integrated_net_GDP_2(closest_index_15,2)

disp('2C, GWP w.r.t. no mitigation, through 2335 3% discount')
integrated_net_GDP_2(closest_index_2,2)

disp('3C, GWP w.r.t. no mitigation, through 2335 3% discount')
integrated_net_GDP_2(closest_index_3,2)

disp('1.5 - 2.0')
integrated_net_GDP_2(closest_index_15,2) - integrated_net_GDP_2(closest_index_2,2)

% int_gdp_min = min(horzcat(min(min(integrated_net_GDP_1)),min(min(integrated_net_GDP_2))));
% int_gdp_max = max(horzcat(max(max(integrated_net_GDP_1)),max(max(integrated_net_GDP_2))));

%% plot

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',7);
set(0,'defaultAxesFontName', 'helvetica')
%colormap(redblue)
%colormap(redgreencmap)
%colormap(pink)

% rgb = [ ...
%    158     1    66  
%    213    62    79   
%    244   109    67   
%    253   174    97   
%    254   224   139   
%    255   255   255   
%    230   245   152   
%    171   221   164
%    102   194   165   
%     50   136   189   
%    94    79   162   ] / 255;

rgb = [ ...
   158     1    66  
   213    62    79   
   244   109    67   
   253   174    97   
   254   224   139 
   255   255   255  
   230   245   152      
   230   245   152   
   230   245   152
   230   245   152         
   230   245   152   
   230   245   152      
   230   245   152 
   171   221   164            
   171   221   164         
   171   221   164      
   171   221   164   
   171   221   164
   171   221   164 
   171   221   164
   102   194   165            
   102   194   165         
   102   194   165      
   102   194   165   
   102   194   165
   102   194   165   
   102   194   165
    50   136   189               
    50   136   189            
    50   136   189         
    50   136   189      
    50   136   189   
    50   136   189
    50   136   189 
   94    79   162                   
   94    79   162               
   94    79   162           
   94    79   162       
   94    79   162   
   94    79   162      
   94    79   162] / 255;

num_contours = 30;
%red_vector = linspace(1,0,num_contours);
%green_vector = linspace(0,1,num_contours);
%blue_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,0,num_contours/2));
%blue_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,0,num_contours/2));

red_vector = horzcat(linspace(1,1,num_contours/2),linspace(1,0,num_contours/2));
green_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,1,num_contours/2));
blue_vector = horzcat(linspace(0,1,num_contours/2),linspace(1,0,num_contours/2));

% figure
% hold on
% plot(1:num_contours,red_vector,'-r')
% plot(1:num_contours,green_vector,'-g')
% plot(1:num_contours,blue_vector,'-b')

rgb2 = horzcat(red_vector',green_vector',blue_vector');

b = repmat(linspace(0,1,200),20,1);
imshow(b,[],'InitialMagnification','fit')
colormap(rgb2)

        subplot(2,2,2)
        hold on
        %grid on
        
        %plot for legend
        
               for discount_rate_i = 1:length(discount_rates)

                    plot(integrated_net_GDP_1(:,discount_rate_i),...
                        max_T(start_temp_ind:end),...
                        discount_rates_styles_1{discount_rate_i},'LineWidth',2)

               end
               
               for discount_rate_i = 1:length(discount_rates)
               
                    plot(integrated_net_GDP_2(:,discount_rate_i),...
                        max_T(start_temp_ind:end),...
                        discount_rates_styles_2{discount_rate_i},'LineWidth',2)
                    
               end
        
              % legend(legend_text_full,'location','NorthEast')
               %legend(legend_text_full,'location','NorthWest')
               
        %plot again with ID of 'best' temperature 
                       
        for discount_rate_i = 1:length(discount_rates)
            
            plot(integrated_net_GDP_1(:,discount_rate_i),...
                max_T(start_temp_ind:end),...
                discount_rates_styles_1{discount_rate_i},'LineWidth',2)
            
            plot(integrated_net_GDP_2(:,discount_rate_i),...
                max_T(start_temp_ind:end),...
                discount_rates_styles_2{discount_rate_i},'LineWidth',2)
            
        end
        
        for discount_rate_i = 1:length(discount_rates)
            
            opt_temp_i_1 = find(squeeze(integrated_net_GDP_1(:,discount_rate_i)) == max(squeeze(integrated_net_GDP_1(:,discount_rate_i))));
            opt_temp_i_2 = find(squeeze(integrated_net_GDP_2(:,discount_rate_i)) == max(squeeze(integrated_net_GDP_2(:,discount_rate_i))));
            
            scatter(integrated_net_GDP_1(opt_temp_i_1,discount_rate_i),...
                    max_T(opt_temp_i_1),...
                    50,...
                    discount_rates_colors{discount_rate_i},...
                    discount_rates_symbols_1{discount_rate_i},...
                    'filled')
                
            scatter(integrated_net_GDP_2(opt_temp_i_2,discount_rate_i),...
                    max_T(opt_temp_i_2),...
                    50,...
                    discount_rates_colors{discount_rate_i},...
                    discount_rates_symbols_2{discount_rate_i},...
                    'filled')

        end
        
        %plot lines
        plot([0 0],[min(max_T) max(max_T)],'-k')
        
%         plot([min([min(min(integrated_net_GDP_1)) min(min(integrated_net_GDP_2))]) max([max(max(integrated_net_GDP_1)) max(max(integrated_net_GDP_2))])],[3 3],'--k')
%         plot([min([min(min(integrated_net_GDP_1)) min(min(integrated_net_GDP_2))]) max([max(max(integrated_net_GDP_1)) max(max(integrated_net_GDP_2))])],[2 2],'--k')
%         plot([min([min(min(integrated_net_GDP_1)) min(min(integrated_net_GDP_2))]) max([max(max(integrated_net_GDP_1)) max(max(integrated_net_GDP_2))])],[1.5 1.5],'--k')
%          plot([-2 6],[3 3],'--k')
%          plot([-2 6],[2 2],'--k')
%          plot([-2 6],[1.5 1.5],'--k')

         plot([-2 6],[3 3],'-m')
         plot([-2 6],[2 2],'-g')
         plot([-2 6],[1.5 1.5],'-b')  
         
        xlabel('net per-capita GWP gain w.r.t. no mitigation')
        ylabel('level of mitigation effort (global warming, C)')
        title('Gross World Product per person, net present value')
        
        ylim([1.4 4])
        %xlim([1.1*int_gdp_min 1.1*int_gdp_max])
        xlim([-2 6])
        
        subplot(2,2,1)
        hold on

                plot_map_vector = reshape(squeeze(two_d_plot(start_time_ind:end,start_temp_ind:end)),[length(time_desc(start_time_ind:end))*length(max_T(start_temp_ind:end)),1]);
                max_color_val = prctile(abs(plot_map_vector),95);

                %max_color_range = 0.8*max(max(abs(two_d_plot_trimmed_2)));
                %min_color_range = -1*max_color_range;
%                 max_color_range = 19;
%                 min_color_range = -3;
                max_color_range = 10;
                min_color_range = -10;

               % num_contours = 30;

                contourf(time_desc(start_time_ind:end)',max_T(start_temp_ind:end),real(two_d_plot(start_time_ind:end,start_temp_ind:end)'),linspace(min_color_range,max_color_range,num_contours), 'LineStyle','none');
                contour(time_desc(start_time_ind:end)',max_T(start_temp_ind:end),real(two_d_plot(start_time_ind:end,start_temp_ind:end)'),[0 0],'k','LineStyle','-','LineWidth',1);
                caxis([min_color_range max_color_range])
                
                plot([end_year_calc_1 end_year_calc_1],[min(max_T) max(max_T)],':k')
                plot([end_year_calc_2 end_year_calc_2],[min(max_T) max(max_T)],':k')
                
%                 plot([2020 end_year_plot],[1.5 1.5],'--k')
%                 plot([2020 end_year_plot],[2 2],'--k')
%                 plot([2020 end_year_plot],[3 3],'--k')
                
                 plot([2020 end_year_plot],[3 3],'-m')
                 plot([2020 end_year_plot],[2 2],'-g')
                 plot([2020 end_year_plot],[1.5 1.5],'-b')  

                %sanePColor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');

                %pcolor(time_desc(3:end)',max_T,two_d_plot(3:end,:)');

                h = colorbar;

                xlim([2021 end_year_plot])
                ylim([1.3 3.75])
                
                %set(h,'YTick',linspace(min_color_range,max_color_range,12))
                
%                 h.Ruler.Scale = 'log';
%                 h.Ruler.MinorTick = 'on';

               % set(gca,'colorscale','log')
               
               
               
               %plot optimal lines
               
                   discount_rate_i = 2;

                   opt_temp_i_1 = find(squeeze(integrated_net_GDP_1(:,discount_rate_i)) == max(squeeze(integrated_net_GDP_1(:,discount_rate_i))));
                   opt_temp_i_2 = find(squeeze(integrated_net_GDP_2(:,discount_rate_i)) == max(squeeze(integrated_net_GDP_2(:,discount_rate_i))));

                   plot([2020 2100],[max_T(opt_temp_i_1) max_T(opt_temp_i_1)],'-c')
                   plot([2020 2160],[max_T(opt_temp_i_2) max_T(opt_temp_i_2)],'-c')


                ylabel(h,ylabels{var_to_plot_i});
                
                ylabel('level of mitigation effort (global warming,C)')
                xlabel('year')
                
                title('Gross World Product per person per year')
                
                
                
%         subplot(2,2,4)
%         hold on
%         %grid on
%         
%             for discount_rate_i = 1:length(discount_rates)
%                 
%                 plot(time_desc(start_time_ind:end_year_calc_ind)',...
%                     squeeze(R(discount_rate_i,1:size(two_d_plot_trimmed_discounted,1))'),...
%                     discount_rates_styles{discount_rate_i},'LineWidth',2)
% 
%             end
%         
% %             plot(time_desc(start_time_ind:end_year_plot_ind)',ones(length(time_desc(start_time_ind:end_year_plot_ind)),1),'-k','LineWidth',2)
% %             plot(time_desc(start_time_ind:end_year_plot_ind)',squeeze(overall_plot_series_all(8,start_time_ind:end_year_plot_ind,1)),'-r','LineWidth',2)
% 
%     %         plot(squeeze(nansum(two_d_plot(start_time_ind:end_year_plot_ind,start_temp_ind:end).*squeeze(overall_plot_series_all(8,start_time_ind:end_year_plot_ind,start_temp_ind:end)),1)),...
%     %             max_T(start_temp_ind:end),'-r','LineWidth',2)
%     %         
%     %         legend('no discount','discount = 3%')
%     %         
%     %         xlabel('integrated GDP loss')
%     %         ylabel('eventual temperature stabilization')
%     %         
%              ylim([0 1.05])
%              xlim([2020 end_year_plot+20])                
                
print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/two_d_per_cap_cons','-dpdf','-bestfit')


%% save stuff for other programs

%pick a time horizon and discount rate

%for 2100, 3% discount rate
    %Nordhaus optimum T = +5K
    %Burke optimim T = +3K

%time_horizon_error_analysis = 2100; this means use integrated_net_GDP_1
%discount_rate_error_analysis = 0.03; this means use the 2nd discount rate index
discount_rate_i = 2;

    % find the aopt index that is optimum
    
    opt_temp_i_1 = find(squeeze(integrated_net_GDP_1(:,discount_rate_i)) == max(squeeze(integrated_net_GDP_1(:,discount_rate_i))));
    aopt_time_series = overall_plot_series_all(14,:,opt_temp_i_1);

code_string = ['_a',num2str(abatement_cost_on),...
                     '_B',num2str(Burke_damage_on),...
                     '_o',num2str(optimize_on),...
                     '_u',num2str(optimize_for_total_utility),...
                     '_uf',num2str(utility_function)];
                 
% strcat('overall_plot_series_',code_string) = overall_plot_series;
% strcat('NPV_utility_',code_string) = NPV_utility;

      save(['/home/pbrown/Documents/MATLAB/DICE/test_Burke/scale_with_temperature_Burke_damage_aopt',code_string,'.mat'],...
          'opt_temp_i_1',...
          'aopt_time_series',...
          'overall_plot_series_all',...
          'R',...
          'titles',...
          'ylabels',...
          'code_string',...
          'end_year_calc_1',...
          'end_year_calc_2',...
          'end_year_plot',...
          'end_year_plot_ind',...
          'end_year_calc_ind_1',...
          'end_year_calc_ind_2',...
          'two_d_plot_trimmed_1',...
          'two_d_plot_trimmed_2',...
          'two_d_plot_trimmed_discounted_1',...
          'two_d_plot_trimmed_discounted_2',...
          'integrated_net_GDP_1',...
          'integrated_net_GDP_2')


%% calculate NPV of mitigation GDP loss to compar with SSP chart

time_2100_i = find(time_desc == 2100);
mean_temp_time_inds = find(time_desc >= 2081 & time_desc <= 2100);

discount_rate_i = 3; %3 is 5%

%load baseline GWP

load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/Run_DICE_no_damage_no_mitigation.mat',...
    'overall_plot_series_all_no_damage_no_mitigation',...
    'titles_no_damage_no_mitigation',...
    'ylabels_no_damage_no_mitigation')

GWP_no_mitigation = overall_plot_series_all_no_damage_no_mitigation(9,:);

GWP_with_mitigation = NaN(1,length(time_desc),num_scens);
mitigation_cost_fraction = overall_plot_series_all(16,:,:);
%mitigation_cost_fraction(1,1,:) = 0;

for exp_i = 1:num_scens
    
    GWP_with_mitigation(1,:,exp_i) = (1-squeeze(mitigation_cost_fraction(1,:,exp_i))).*squeeze(GWP_no_mitigation);
    
end

mitigation_cost = GWP_with_mitigation - GWP_no_mitigation;

mitigation_cost_discounted_5 = NaN(1,length(time_desc),num_scens);
GWP_no_mitigation_discounted_5 = R(discount_rate_i,:).*squeeze(GWP_no_mitigation);

mitigation_cost_frac_NPV = NaN(1,num_scens);
mean_temp_2081_2100 = NaN(1,num_scens);

for exp_i = 1:num_scens
    
    mitigation_cost_discounted_5(1,:,exp_i) = R(discount_rate_i,:).*squeeze(mitigation_cost(1,:,exp_i));
    
    mitigation_cost_frac_NPV(exp_i) = nansum(squeeze(-1*mitigation_cost_discounted_5(1,1:time_2100_i,exp_i)))./nansum(GWP_no_mitigation_discounted_5(1:time_2100_i));
    
    mean_temp_2081_2100(exp_i) = mean(squeeze(overall_plot_series_all(12,mean_temp_time_inds,exp_i)));
    
end

% plot NPV mitigation cost vs. global temperature 

data_xmin = 0;
data_xmax = 1.1*nanmax(mitigation_cost_frac_NPV);
data_ymin = 0.9*nanmin(mean_temp_2081_2100);
data_ymax = 1.1*nanmax(mean_temp_2081_2100);

rcp_26_45_60_85_temps = [1.6 2.4 2.8 4.3];
rcp_26_45_60_85_temps_sig = [0.4 0.5 0.5 0.7];

rcp_26_45_60_85_temps_maxes = rcp_26_45_60_85_temps + rcp_26_45_60_85_temps_sig;
rcp_26_45_60_85_temps_mins = rcp_26_45_60_85_temps - rcp_26_45_60_85_temps_sig;

rcp_26_45_60_85_NVP_mitig_SSPall_maxes = (1/100)*[6.96 1.87 0.69 0.01];
rcp_26_45_60_85_NVP_mitig_SSPall_mins = (1/100)*[0.38 -0.01 0.01 -0.01];

rcp_26_45_60_85_NVP_mitig_SSPall_half_range = (rcp_26_45_60_85_NVP_mitig_SSPall_maxes - rcp_26_45_60_85_NVP_mitig_SSPall_mins)./2;
rcp_26_45_60_85_NVP_mitig_SSPall_central = rcp_26_45_60_85_NVP_mitig_SSPall_half_range+rcp_26_45_60_85_NVP_mitig_SSPall_mins;



rcp_26_45_60_NVP_mitig_SSP5 = (1/100)*[2.34 1.28 0.69 0];
rcp_26_45_60_NVP_mitig_SSP5_sig = (1/100)*[(2.91-2.34)/4 ...
                                   (1.87-0.86)/4 ...
                                   (0.69-0.56)/4 ...
                                   0.1];
%                                
% rcp_26_45_60_NVP_mitig_SSP2 = (1/100)*[0.76 0.14 0.03];
% rcp_26_45_60_NVP_mitig_SSP2_sig = (1/100)*[(2.91-2.34)/4 ...
%                                    (1.87-0.86)/4 ...
%                                    (0.69-0.56)/4];
rcp_colors_lines = {'-g','-b','-m','-r'};
rcp_colors = {'g','b','m','r'};

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',12);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
    for rcp_i = 1:length(rcp_26_45_60_85_temps)
        boundedline([rcp_26_45_60_85_NVP_mitig_SSPall_mins(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_maxes(rcp_i)],...
        [rcp_26_45_60_85_temps(rcp_i) rcp_26_45_60_85_temps(rcp_i)],...
        [rcp_26_45_60_85_temps_sig(rcp_i) rcp_26_45_60_85_temps_sig(rcp_i)],...
        rcp_colors_lines{rcp_i},'alpha','transparency',0.1)
    end
    
    for rcp_i = 1:length(rcp_26_45_60_85_NVP_mitig_SSPall_central)
        boundedline([rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i)],...
        [rcp_26_45_60_85_temps_mins(rcp_i) rcp_26_45_60_85_temps_maxes(rcp_i)],...
        [rcp_26_45_60_85_NVP_mitig_SSPall_half_range(rcp_i) rcp_26_45_60_85_NVP_mitig_SSPall_half_range(rcp_i)],...
        rcp_colors_lines{rcp_i},'orientation', 'horiz','alpha','transparency',0.1)
    end    
    
     plot(mitigation_cost_frac_NPV,mean_temp_2081_2100,'-k','LineWidth',2)
     scatter(mitigation_cost_frac_NPV,mean_temp_2081_2100,25,[0 0 0],'filled')
     
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
     for rcp_i = 1:length(rcp_26_45_60_NVP_mitig_SSP5)
        scatter(rcp_26_45_60_85_NVP_mitig_SSPall_central(rcp_i),rcp_26_45_60_85_temps(rcp_i),100,rcp_colors{rcp_i},'o','filled')
     end
  
         %xlim([2015 end_year_plot])
 
         xlabel('mitigation cost: GWP losses 2010-2100, discounted at 5% (fraction)')
         ylabel('global temp 2081-2100 (C)')
         %title(titles{var_i})    

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/NPV_mitig_vs_avg_temp_2081_2100','-dpdf','-bestfit')

%% calculate an estimate of the raw Burke damage

good_time_ind = find(time_desc == 2100);

%the S that exists here is the no mitigation S

GWP_with_damage_no_mitigation = S(20,2:good_time_ind+1);
GWP_baseline_no_damage = S(24,1:good_time_ind);

disp('2100 DICE baseline GWP:')
GWP_baseline_no_damage(end)

disp('2100 SSP1 baseline GWP')
SSP1_2100_baseline_SSP = 565.389; %https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=series

disp('scaling factor')
scale_GDP_ratio = 565.389./GWP_baseline_no_damage(end);

%scale_GDP_ratio = 1;

%scale dice time series to SSP1 and take difference

diff_time_series = GWP_baseline_no_damage.*scale_GDP_ratio - GWP_with_damage_no_mitigation.*scale_GDP_ratio;

discount_rate_i = 2;

diff_time_series_discount = squeeze(R(discount_rate_i,1:good_time_ind)).*diff_time_series;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',12);
set(0,'defaultAxesFontName', 'helvetica')

    subplot(2,2,1)
    hold on
    
        plot(time_desc(1:good_time_ind),GWP_baseline_no_damage,'-k','LineWidth',2)
        plot(time_desc(1:good_time_ind),GWP_with_damage_no_mitigation,'-r','LineWidth',2)
        
        title('GWP')
        
    subplot(2,2,2)
    hold on
    
        plot(time_desc(1:good_time_ind),diff_time_series,'-k','LineWidth',2)
        
        title('difference')
        
    subplot(2,2,3)
    hold on
    
        plot(time_desc(1:good_time_ind),R(discount_rate_i,1:good_time_ind),'-k','LineWidth',2)
        
        title('discounted rate')        
        
    subplot(2,2,4)
    hold on
    
        plot(time_desc(1:good_time_ind),diff_time_series_discount,'-k','LineWidth',2)
        
        title('discounted difference')
        
disp('Dice-Burke GWP damage (NPV, 2100, 3%)')

diff_time_series_disc = nansum(squeeze(R(discount_rate_i,1:good_time_ind)).*diff_time_series)