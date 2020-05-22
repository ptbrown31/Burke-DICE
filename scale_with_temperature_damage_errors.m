close all
clear all

%% setpath for rest of specific functions

addpath('/home/pbrown/Documents/MATLAB/DICE/test_Burke/DICE_Burke_Damage_growth')

% nordhaus
load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/scale_with_temperature_Burke_damage_aopt_a1_B0_o0_u1_uf1.mat',...
  'opt_temp_i_1',...
  'aopt_time_series',...
  'overall_plot_series_all',...
  'R',...
  'titles',...
  'ylabels',...
  'code_string',...
  'end_year_plot',...
  'end_year_calc_1',...
  'end_year_calc_2',...
  'end_year_plot_ind',...
  'end_year_calc_ind_1',...
  'end_year_calc_ind_2',...
  'two_d_plot_trimmed_1',...
  'two_d_plot_trimmed_2',...
  'two_d_plot_trimmed_discounted_1',...
  'two_d_plot_trimmed_discounted_2',...
  'integrated_net_GDP_1',...
  'integrated_net_GDP_2')

opt_temp_i_1_nordhaus = opt_temp_i_1;
overall_plot_series_all_nordhaus = overall_plot_series_all;

integrated_net_GDP_1_nordhaus = integrated_net_GDP_1;
integrated_net_GDP_2_nordhaus = integrated_net_GDP_2;

integrated_net_GDP_nordhaus_both_time_horz = cat(3,integrated_net_GDP_1_nordhaus,integrated_net_GDP_2_nordhaus);

%burke
load('/home/pbrown/Documents/MATLAB/DICE/test_Burke/scale_with_temperature_Burke_damage_aopt_a1_B1_o0_u1_uf1.mat',...
  'opt_temp_i_1',...
  'aopt_time_series',...
  'overall_plot_series_all',...
  'R',...
  'titles',...
  'ylabels',...
  'code_string',...
  'end_year_plot',...
  'end_year_calc_1',...
  'end_year_calc_2',...
  'end_year_plot_ind',...
  'end_year_calc_ind_1',...
  'end_year_calc_ind_2',...
  'two_d_plot_trimmed_1',...
  'two_d_plot_trimmed_2',...
  'two_d_plot_trimmed_discounted_1',...
  'two_d_plot_trimmed_discounted_2',...
  'integrated_net_GDP_1',...
  'integrated_net_GDP_2')

opt_temp_i_1_burke = opt_temp_i_1;
overall_plot_series_all_burke = overall_plot_series_all;

integrated_net_GDP_1_burke = integrated_net_GDP_1;
integrated_net_GDP_2_burke = integrated_net_GDP_2;

integrated_net_GDP_burke_both_time_horz = cat(3,integrated_net_GDP_1_burke,integrated_net_GDP_2_burke);

integrated_net_GDP_both_both = cat(4,integrated_net_GDP_nordhaus_both_time_horz,integrated_net_GDP_burke_both_time_horz);

% size(integrated_net_GDP_both_both)
% ans =
% 
%     64     3     2     2

% 1) - level of mitigation (assume nordhause or burke)
% 2) - discount rate
% 3) - time horizon
% 4) - get nordhaus or burke

table_assume_nordhaus_get_nordhaus = squeeze(integrated_net_GDP_both_both(opt_temp_i_1_nordhaus,:,:,1));
table_assume_burke_get_burke = squeeze(integrated_net_GDP_both_both(opt_temp_i_1_burke,:,:,2));
table_assume_nordhaus_get_burke = squeeze(integrated_net_GDP_both_both(opt_temp_i_1_nordhaus,:,:,2));
table_assume_burke_get_nordhaus = squeeze(integrated_net_GDP_both_both(opt_temp_i_1_burke,:,:,1));


%% make combo arrays

all_assume_nordhaus_get_nordhaus = squeeze(overall_plot_series_all_nordhaus(:,:,opt_temp_i_1_nordhaus));
all_assume_burke_get_burke = squeeze(overall_plot_series_all_burke(:,:,opt_temp_i_1_burke));
all_assume_nordhaus_get_burke = squeeze(overall_plot_series_all_burke(:,:,opt_temp_i_1_nordhaus));
all_assume_burke_get_nordhaus = squeeze(overall_plot_series_all_nordhaus(:,:,opt_temp_i_1_burke));

all_all = cat(3,all_assume_nordhaus_get_nordhaus,...
              all_assume_burke_get_burke,...
              all_assume_nordhaus_get_burke,...
              all_assume_burke_get_nordhaus);

%% plot tight

linstyles = {'-g','-.r',':k','--b'};

time_desc = 2010:5:2330;

sub_var_inds = [10 12 16 13 2 4];

row_num = 3;
col_num = 2;

end_year_plot = 2160;

    FigHandle = figure('Position', [100, 100, 1100, 600]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',6);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    grid on
    
for var_i = 1:length(sub_var_inds)
    
    subplot(row_num,col_num,var_i)
    hold on
    %grid on
    
        for exp_i = 1:4
    
            if var_i ~= 4; plot(time_desc,squeeze(all_all(sub_var_inds(var_i),:,exp_i)),linstyles{exp_i},'LineWidth',2); end
            if var_i == 4; plot(time_desc,-1*squeeze(all_all(sub_var_inds(var_i),:,exp_i)-1),linstyles{exp_i},'LineWidth',2); end
        
        end
        
%         if var_i == 2 || var_i == 3
%             
%             ylim([0 0.25])
%             %set(gca, 'YScale', 'log')
%             
%         end
        
        if var_i == 4
            
            legend('assume Nordhaus, get Nordhaus',...
                   'assume Burke, get Burke',...
                   'assume Nordhaus, get Burke',...
                   'assume Burke, get Nordhaus',...
                   'Location','NorthWest')
            
        end
              
        xlim([2010 end_year_plot])    

        %xlim([2015 end_year_plot])

        xlabel('year')
        ylabel(ylabels{sub_var_inds(var_i)})
        title(titles{sub_var_inds(var_i)})

end

print(gcf,'/home/pbrown/Documents/MATLAB/DICE/test_Burke/guess_get_tight_time_series','-dpdf')


