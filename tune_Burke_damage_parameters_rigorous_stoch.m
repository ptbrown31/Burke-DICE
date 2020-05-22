close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/')

%% read in Burke data
row_offset = 0;
col_offset = 0;

SSPs = {'SSP1','SSP2','SSP3','SSP4','SSP5'};
rcps = {'RCP26','RCP45','RCP60','RCP85'};
time_horzs = {'2049','2099'};

filename_base = '/home/pbrown/Documents/MATLAB/DICE/test_Burke/BDD18_Fig_4a_data/BDD18_Fig4a_GDP_per_cap_vs_GMSAT_';

GDP_loss_Burke_Fig_4a_all_data = NaN(39,2,4,5,2); %1) number of points (models), 2) GMSAT vs. per capita GDP loss, 3)RCPs, 4) SSPs, 5) time integration

for SSPs_i = 1:length(SSPs)
    for rcp_i = 1:length(rcps)
         for time_horizs_i = 1:length(time_horzs)

            filename = strcat(filename_base,SSPs{SSPs_i},'_',rcps{rcp_i},'_',time_horzs{time_horizs_i},'.csv');
            
            d = dir(filename);
            
            if isempty(d) == 0

                M = csvread(filename,row_offset,col_offset);

                num_rows = size(M,1);

                GDP_loss_Burke_Fig_4a_all_data(1:num_rows,:,rcp_i,SSPs_i,time_horizs_i) = M;
                
            end
        end
    end
end

zero_inds = find(GDP_loss_Burke_Fig_4a_all_data == 0);

GDP_loss_Burke_Fig_4a_all_data(zero_inds) = NaN;

%% read in SSP data 

init_time = horzcat(2005,2010:10:2100);

%IAM: IMAGE

RCPs = {'RCP-26',...
       'RCP-34',...
       'RCP-45',...
       'RCP-60',...
       'RCP-Baseline'};
            
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

SSP_all_data_raw = NaN(length(RCPs),length(init_time),length(var_names),length(SSPs)); %1) SSP1 sub scenarios (RCPs), 2) time, 3) variable

for SSPs_i = 1:length(SSPs)
    
    %read
    SSP_all_data_raw(:,:,1,SSPs_i) = xlsread(filename,SSPs{SSPs_i},'E5:O9');
    SSP_all_data_raw(:,:,2,SSPs_i) = xlsread(filename,SSPs{SSPs_i},'E11:O15');
    SSP_all_data_raw(:,:,3,SSPs_i) = xlsread(filename,SSPs{SSPs_i},'E17:O21');
    SSP_all_data_raw(:,:,4,SSPs_i) = xlsread(filename,SSPs{SSPs_i},'E23:O27');
    SSP_all_data_raw(:,:,5,SSPs_i) = xlsread(filename,SSPs{SSPs_i},'E29:O33');
    
end

no_data_i = find(SSP_all_data_raw == -99999);
SSP_all_data_raw(no_data_i) = NaN;

%interpolate in time

time_desc = 2005:5:2100;

time_desc_2050_i = find(time_desc == 2050);
time_desc_2100_i = find(time_desc == 2100);

SSP_all_data = NaN(length(RCPs),length(time_desc),length(var_names),length(SSPs));
for SSPs_i = 1:length(SSPs)
    for rcp_i = 1:length(RCPs)
        for var_names_i = 1:length(var_names)

            SSP_all_data(rcp_i,:,var_names_i,SSPs_i) = interp1(init_time,...
                                                      squeeze(SSP_all_data_raw(rcp_i,:,var_names_i,SSPs_i)),...
                                                      time_desc);
        end
    end
end

%make a gdp per capita 

gdp_per_capita = SSP_all_data(:,:,1,:)./SSP_all_data(:,:,2,:);

SSP_all_data = cat(3,SSP_all_data,gdp_per_capita);

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
         
for SSPs_i = 1:length(SSPs)
    for rcp_i = 1:length(RCPs)
        for time_i = 1:length(time_desc)

            SSP_all_data(rcp_i,time_i,:,SSPs_i) = conversion_farctors'.*squeeze(SSP_all_data(rcp_i,time_i,:,SSPs_i));

        end
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
    A0 = 5.115; %3.8; % DICE-2013R; 27.22 for DICE-2007
    %   Initial growth rate for TFP per 5 years
    Ag0 = 0.076; %0.079; % DICE-2013R
    %   Decline rate of TFP per 5 years
    Agd = 0.005; %0.006; % DICE-2013R
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
    
%% main loop

num_values_test = 100;

        %choose a central coeff
        damage_TFP_growth_coeffs_on_T = 0.01*rand(num_values_test,1);
        damage_dep_rate_coeffs_on_T = 0.1*rand(num_values_test,1);
        
%         damage_TFP_growth_coeffs_on_T = range_values_1.*0.0002;
%         damage_dep_rate_coeffs_on_T = range_values_2.*0.13;        
%             damage_TFP_growth_coeffs_on_T = range_values*0.0024;
%             damage_dep_rate_coeffs_on_T = range_values*0.04;
            
errors_between_my_calc_and_Burke_fit = NaN(length(RCPs),length(SSPs),2,num_values_test,num_values_test);
rmse_across_all_SSPs = NaN(num_values_test,num_values_test);

%main results (black squares)
Y_per_Ls_vs_no_damage = NaN(length(RCPs),length(time_desc),length(SSPs),num_values_test,num_values_test);
Y_per_Ls_vs_no_damage_points = NaN(length(RCPs),length(time_horzs),length(SSPs),num_values_test,num_values_test);

lin_count = 1;

for dep_rate_value_i = 1:num_values_test
    for TFP_growth_rate_value_i = 1:num_values_test

%     dep_rate_value_i = 10;
%     TFP_growth_rate_value_i = 10;
        
        lin_count./(num_values_test.^2);
        
        %first select the coefficient increasing capital depreciation and
        %reducing TFP growth

        damage_dep_rate_coeff_on_T = damage_dep_rate_coeffs_on_T(dep_rate_value_i);
        damage_TFP_growth_coeff_on_T = damage_TFP_growth_coeffs_on_T(TFP_growth_rate_value_i);
        
        
        %define empties for this paramater combo
        
            %I have 4 global temperature trajectories (RCPs)
            %this are connected to TFP and K depreciation
            
            A_with_damage = NaN(length(RCPs),length(time_desc),length(SSPs));
            K_with_damage = NaN(length(RCPs),length(time_desc),length(SSPs));
            Y_with_damage = NaN(length(RCPs),length(time_desc),length(SSPs));
            Y_per_L_with_damage = NaN(length(RCPs),length(time_desc),length(SSPs));

            Y_without_damage = NaN(length(RCPs),length(time_desc),length(SSPs));
            Y_per_L_without_damage = NaN(length(RCPs),length(time_desc),length(SSPs));            

            %need to run model forward in time such that it produces correct losses

            deltak = 0.1;
            %optlrsav = 0.258278145695364;



        %make damage calculations for different SSP-RCP combinations
        for SSPs_i = 1:length(SSPs)
            for rcp_i = 1:length(RCPs)    

                % have TFP background (from DICE), population (from SSP), gross output
                % (from SSP), need to back out kapital that is consisten with these for
                % this given SSP

                %rcp_i = 1;

                %http://www.wolframalpha.com/input/?i=solve+for+K,+Y%3DA*(K%5Eg)*(L%5E(1-g))

                Y = squeeze(SSP_all_data(rcp_i,:,1,SSPs_i));
                L = squeeze(SSP_all_data(rcp_i,:,2,SSPs_i));

                %only continue if there is data for this SSP-RCP combo
                
                if sum(isnan(Y)) == 0

                    % this is what the production function looks like: S2(7) = A_damaged * (S2(1) ^ gama) * (L(t1 + 1)/1000)^(1 - gama);

                    K = ((Y.*((L./1000).^(gama-1)))./(A)).^(1./gama);

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
                        
%                         figure
%                         plot(save_rate)
                        
            %% calculate damages, accumulating over time
                
                    A_with_damage(rcp_i,1,SSPs_i) = A(1); %same initial value goes in for each RCP
                    K_with_damage(rcp_i,1,SSPs_i) = K(1);
                    Y_with_damage(rcp_i,1,SSPs_i) = Y(1);
                    Y_per_L_with_damage(rcp_i,1,SSPs_i) = Y(1)./squeeze(SSP_all_data_raw(rcp_i,1,2,SSPs_i));

                    Y_per_L_without_damage(rcp_i,1,SSPs_i) = Y(1)./squeeze(SSP_all_data_raw(rcp_i,1,2,SSPs_i));
                    Y_without_damage(rcp_i,1,SSPs_i) = Y(1)./squeeze(SSP_all_data_raw(rcp_i,1,2,SSPs_i));

                        for time_i = 2:length(time_desc)

                            %what is the global warming?

                            global_temp_now = SSP_all_data(rcp_i,time_i,3,SSPs_i);

                            %calculate the damaged TFP for this temperature

                    %                Ag_damaged = Ag(t1-1) - damage_TFP_growth_coeff_on_T*S1(2);
                    %                A_damaged = S1(23)/(1-Ag_damaged);  %S1(23) is the previous timesteps damaged TFP

                                     Ag_damaged_now = Ag(time_i-1) - damage_TFP_growth_coeff_on_T*global_temp_now;

                                     A_with_damage(rcp_i,time_i,SSPs_i) = A_with_damage(rcp_i,time_i-1,SSPs_i)./(1-Ag_damaged_now);
                    %                 A_with_damage(ssp_sub_scen_i,time_i) = A(time_i);

                            %calculate the damaged deltak for this temperature

                    %               deltak_damaged = damage_dep_rate_coeff_on_T*S1(2);
                    %               if deltak_damaged < deltak; deltak_damaged = deltak; end
                    
                                    deltak_damaged_now = damage_dep_rate_coeff_on_T*global_temp_now;
                                    if deltak_damaged_now < deltak; deltak_damaged_now = deltak; end

                                    %deltak_damaged_now = deltak + damage_dep_rate_coeff_on_T*global_temp_now;
                                    %if deltak_damaged_now < deltak; deltak_damaged_now = deltak; end

                            %calculate the damaged capital stock (dont use optlrsave but
                            %emperically derived time-varying saveings rate

                    %               S2(1) = tstep * optlrsav * Q + (1 - deltak_damaged) ^ tstep * S1(1);

                    %                 K_with_damage(ssp_sub_scen_i,time_i) = tstep * optlrsav * Y_with_damage(ssp_sub_scen_i,time_i-1) + ...
                    %                                                        (1 - deltak_damaged_now) ^ tstep * K_with_damage(ssp_sub_scen_i,time_i-1);

                                    K_with_damage(rcp_i,time_i,SSPs_i) = tstep * save_rate(time_i-1) * Y_with_damage(rcp_i,time_i-1,SSPs_i) + ...
                                                                           (1 - deltak_damaged_now) ^ tstep * K_with_damage(rcp_i,time_i-1,SSPs_i);

                            %calculate the damaged gross output

                    %               S2(7) = A_damaged * (S2(1) ^ gama) * (L(t1 + 1)/1000)^(1 - gama);

                                    Y_with_damage(rcp_i,time_i,SSPs_i) = A_with_damage(rcp_i,time_i,SSPs_i) * ...
                                                                          (K_with_damage(rcp_i,time_i,SSPs_i) ^ gama) * ...
                                                                          (SSP_all_data(rcp_i,time_i,2,SSPs_i)/1000)^(1 - gama);

                            %calculate undamaged gross output

                                    Y_without_damage(rcp_i,time_i,SSPs_i) = A(time_i) * ...
                                                                              (K(time_i) ^ gama) * ...
                                                                              (SSP_all_data(rcp_i,time_i,2,SSPs_i)/1000)^(1 - gama);

                            %calculate the damaged gross output per capita

                                    Y_per_L_with_damage(rcp_i,time_i,SSPs_i) = Y_with_damage(rcp_i,time_i,SSPs_i)./SSP_all_data(rcp_i,time_i,2,SSPs_i);

                                    Y_per_L_without_damage(rcp_i,time_i,SSPs_i) = Y_without_damage(rcp_i,time_i,SSPs_i)./SSP_all_data(rcp_i,time_i,2,SSPs_i);

                                    Y_per_Ls_vs_no_damage(rcp_i,time_i,SSPs_i,dep_rate_value_i,TFP_growth_rate_value_i) = -100*(1-(Y_per_L_with_damage(rcp_i,time_i,SSPs_i)./Y_per_L_without_damage(rcp_i,time_i,SSPs_i)));

                        end

            Y_per_Ls_vs_no_damage_points(:,:,:,:,:) = squeeze(Y_per_Ls_vs_no_damage(:,[time_desc_2050_i time_desc_2100_i],:,:,:));
            Temp_points = squeeze(SSP_all_data(:,[time_desc_2050_i time_desc_2100_i],3,:));
            
            
            %fit function to the Burke data
            
                    for time_horizs_i = 1:length(time_horzs)

                        x = reshape(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,1,:,SSPs_i,time_horizs_i)),...
                            numel(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,1,:,SSPs_i,time_horizs_i))),1);

                        y = reshape(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,2,:,SSPs_i,time_horizs_i)),...
                            numel(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,2,:,SSPs_i,time_horizs_i))),1);

                        x = x(isnan(x) == 0);
                        y = y(isnan(y) == 0);

                            if isempty(x) == 0 && isempty(y) == 0

                                 myfittype = fittype('a + b*log(x)',...
                                'dependent',{'y'},'independent',{'x'},...
                                'coefficients',{'a','b'});

                                myfit = fit(x,y,myfittype);


                                        Burke_damage_at_this_Temp = myfit.a + myfit.b*log(Temp_points(rcp_i,time_horizs_i,SSPs_i));

                                        my_damage_at_this_Temp = Y_per_Ls_vs_no_damage_points(rcp_i,time_horizs_i,SSPs_i,dep_rate_value_i,TFP_growth_rate_value_i);

                                        errors_between_my_calc_and_Burke_fit(rcp_i,SSPs_i,time_horizs_i,dep_rate_value_i,TFP_growth_rate_value_i) = ...
                                            (my_damage_at_this_Temp - Burke_damage_at_this_Temp).^2;

                            end
                    end
                end
            end
        end
        
        rmse_across_all_SSPs(dep_rate_value_i,TFP_growth_rate_value_i) = ...
            sqrt(nanmean(nanmean(nanmean(squeeze(errors_between_my_calc_and_Burke_fit(:,:,:,dep_rate_value_i,TFP_growth_rate_value_i))))));
        
    end
end

%% find index of lowest RMSE and plot results

[row_i,col_i] = find(rmse_across_all_SSPs == min(min(min(rmse_across_all_SSPs))));

best_damage_dep_rate_coeffs_on_T = damage_dep_rate_coeffs_on_T(row_i)
best_damage_TFP_growth_coeffs_on_T = damage_TFP_growth_coeffs_on_T(col_i)

save('/home/pbrown/Documents/MATLAB/DICE/test_Burke/tune_Burke_damage_parameters_rigorous_stoch.mat',...
          'GDP_loss_Burke_Fig_4a_all_data',...
          'Temp_points',...
          'Y_per_Ls_vs_no_damage_points',...
          'damage_dep_rate_coeffs_on_T',...
          'damage_TFP_growth_coeffs_on_T',...
          'rmse_across_all_SSPs',...
          'row_i',...
          'col_i',...
          'best_damage_dep_rate_coeffs_on_T',...
          'best_damage_TFP_growth_coeffs_on_T',...
          'SSPs',...
          'time_horzs',...
          'rcps')

        %% plot
        
%         row_num = 3;
%         col_num = 2;
%         
%         colors = {'g','r','b','m'};
%         syms = {'o','o'};
%         
%         FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
%         set(gcf,'color',[1 1 1]);
%         set(0, 'DefaultAxesFontSize',9);
%         set(0,'defaultAxesFontName', 'helvetica')
%         
%         for SSPs_i = 1:length(SSPs)
%             subplot(row_num,col_num,SSPs_i)
%             hold on
%             
%                 for time_horizs_i = 1:length(time_horzs)
%                     for rcp_i = 1:length(rcps)
%                             x = GDP_loss_Burke_Fig_4a_all_data(:,1,rcp_i,SSPs_i,time_horizs_i);
%                             y = GDP_loss_Burke_Fig_4a_all_data(:,2,rcp_i,SSPs_i,time_horizs_i);
% 
%                             x = x(isnan(x) == 0);
%                             y = y(isnan(y) == 0);
% 
%                         if isempty(x) == 0 && isempty(y) == 0
% 
%                             scatter(x,y,50,colors{rcp_i},syms{time_horizs_i},'filled');
% 
%                         end
%                     end
%                 end
%                 
%                     %fit log over all values
%                     
%                 for time_horizs_i = 1:length(time_horzs)
% 
%                         x = reshape(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,1,:,SSPs_i,time_horizs_i)),...
%                             numel(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,1,:,SSPs_i,time_horizs_i))),1);
% 
%                         y = reshape(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,2,:,SSPs_i,time_horizs_i)),...
%                             numel(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,2,:,SSPs_i,time_horizs_i))),1);
% 
%                         x = x(isnan(x) == 0);
%                         y = y(isnan(y) == 0);
% 
%                         if isempty(x) == 0 && isempty(y) == 0
% 
%                              myfittype = fittype('a + b*log(x)',...
%                             'dependent',{'y'},'independent',{'x'},...
%                             'coefficients',{'a','b'});
% 
%                             myfit = fit(x,y,myfittype);
% 
%                             plot(myfit)
% 
%                          scatter(Temp_points(:,time_horizs_i,SSPs_i),...
%                                  Y_per_Ls_vs_no_damage_points(:,time_horizs_i,SSPs_i,row_i,col_i),...
%                                  100,'k','s','filled')
%                              
%                          legend('off')
%                              
%                         end
%                 end
%                     
%                     
%                     rmse = sqrt(nanmean(nanmean(squeeze(errors_between_my_calc_and_Burke_fit(:,SSPs_i,:,row_i,col_i)))));
% 
%                     ylabel('GWP loss (%)')
%                     xlabel('Global temperature above preindustrial (C)')
%                     title(strcat(SSPs{SSPs_i},', RMSE=',num2str(round(rmse,2))))
%                     
%         end
       
         %% plot
%         
%         SSPs_i = 1;
%         
%         colors = {'g','r','b','m'};
%         syms = {'o','o'};
%         
%         FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
%         set(gcf,'color',[1 1 1]);
%         set(0, 'DefaultAxesFontSize',12);
%         set(0,'defaultAxesFontName', 'helvetica')
%         hold on
% 
%                 for time_horizs_i = 1:length(time_horzs)
%                     for rcp_i = 1:length(rcps)
%                             x = GDP_loss_Burke_Fig_4a_all_data(:,1,rcp_i,SSPs_i,time_horizs_i);
%                             y = GDP_loss_Burke_Fig_4a_all_data(:,2,rcp_i,SSPs_i,time_horizs_i);
% 
%                             x = x(isnan(x) == 0);
%                             y = y(isnan(y) == 0);
% 
%                         if isempty(x) == 0 && isempty(y) == 0
%                             
%                             scatter(x,y,50,...
%                             'MarkerFaceColor',colors{rcp_i},...
%                             'MarkerEdgeColor','none',...
%                             'MarkerFaceAlpha',0.4)
% 
%                         end
%                     end
%                 end
%                 
%                     %fit log over all values
%                     
%                 for time_horizs_i = 1:length(time_horzs)
% 
%                         x = reshape(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,1,:,SSPs_i,time_horizs_i)),...
%                             numel(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,1,:,SSPs_i,time_horizs_i))),1);
% 
%                         y = reshape(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,2,:,SSPs_i,time_horizs_i)),...
%                             numel(squeeze(GDP_loss_Burke_Fig_4a_all_data(:,2,:,SSPs_i,time_horizs_i))),1);
% 
%                         x = x(isnan(x) == 0);
%                         y = y(isnan(y) == 0);
% 
%                         if isempty(x) == 0 && isempty(y) == 0
% 
%                              myfittype = fittype('a + b*log(x)',...
%                             'dependent',{'y'},'independent',{'x'},...
%                             'coefficients',{'a','b'});
% 
%                             myfit = fit(x,y,myfittype)
% 
%                             %plot(myfit,x,y)
%                             plot(myfit)
% 
%                          scatter(Temp_points(:,time_horizs_i,SSPs_i),...
%                                  Y_per_Ls_vs_no_damage_points(:,time_horizs_i,SSPs_i,row_i,col_i),...
%                                  100,'k','s','filled')
%                              
%                          legend('off')
%                              
%                         end
%                 end
%                     
%                     
%                     rmse = sqrt(nanmean(nanmean(squeeze(errors_between_my_calc_and_Burke_fit(:,SSPs_i,:,row_i,col_i)))));
% 
%                     ylabel('Gross World Product loss (%)')
%                     xlabel('Global temperature above preindustrial (C)')
%                     title(strcat(SSPs{SSPs_i},', RMSE=',num2str(round(rmse,2))))
%                     
%                     ylim([-40 0])
%                     
%                     print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/otune_Burke_damage_rigorous','-dpdf','-bestfit')



