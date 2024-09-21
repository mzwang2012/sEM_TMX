clear;
close all;

choose_post = 1; % 1 for plot coeff v.s. t
                 % 2 for plot coeff v.s. global mean T

%% load coefficients
load('pca_historical/modes_Tmx_365skip5_global_historical.mat','lambda');
filename = strcat('pca_historical/Tmx_365_global_historical_coeff_',num2str(1,'%02d'),'.mat');
load(filename);
coeff_historical = coeff;

filename = strcat('pca_historical/Tmx_365_global_ssp126_coeff_',num2str(1,'%02d'),'.mat');
load(filename);
coeff_ssp126 = coeff;
filename = strcat('pca_historical/Tmx_365_global_ssp585_coeff_',num2str(1,'%02d'),'.mat');
load(filename);
coeff_ssp585 = coeff;

% load('../ERA5/PCA/Tmx_365_global_1979_2019_coeff.mat');
% coeff_era5 = coeff;
% clear coeff;

n_days_historical = (2015 - 1850)*365;
time_historical = 1:n_days_historical;
n_days_ssp = (2100-2015+1)*365;
time_ssp = (n_days_historical+1):(n_days_historical+n_days_ssp);

lambda = lambda / (n_days_historical/5);
%% plot coeff v.s. t
if(choose_post == 1)
    for i_mode = 5%:20

    ticks = 1:365:(n_days_historical+n_days_ssp);
    ticklabel = num2cell(1850:20:2100);
    for i=1:length(ticklabel)
        ticklabel{i} = num2str(ticklabel{i});
    end

    figure;
    signal = movmean(coeff_historical(:,i_mode)/sqrt(lambda(i_mode)),365);
%     signal = coeff_historical(:,i_mode)/sqrt(lambda(i_mode));
    plot(time_historical,signal,'k','linewidth',1.2);
    hold on;
    signal = movmean(coeff_ssp126(:,i_mode)/sqrt(lambda(i_mode)),365);
%     signal = coeff_ssp126(:,i_mode)/sqrt(lambda(i_mode));
    plot(time_ssp,signal,'g','linewidth',1.2);
    signal = movmean(coeff_ssp585(:,i_mode)/sqrt(lambda(i_mode)),365);
%     signal = coeff_ssp585(:,i_mode)/sqrt(lambda(i_mode));
    plot(time_ssp,signal,'r','linewidth',1.2);

    % istart = (1979-1850)*365+1;
    % iend = (2019 - 1850 +1)*365;
    % signal_era5 = movmean(coeff_era5(:,i_mode)/sqrt(lambda(i_mode)),365);
    % plot(istart:iend,signal_era5,'Color',[0.5 0.5 0.5],'linewidth',1.2);

    xlim([1 n_days_ssp+n_days_historical]);
    % ylim([-4 4]);
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'XTick',ticks(1:20:end),...
        'XTickLabel',ticklabel,'fontsize',15);
    xlabel('time','fontsize',15);
    ylabel('$\bar{a}_i / \sqrt{\lambda_i}$','interpreter','latex','fontsize',20);
    set(gcf,'Units','inches',...
        'Position',[1 1 12 3], ...
        'PaperPosition',[0 0 12.5 3.5]);
    legend('historical','ssp126','ssp585','location','northwest');
%     legend('fontsize',20);
    legend('boxoff');

%     outname = strcat('1011_pca_hist_coeff_',num2str(i_mode,'%3d'));
%     print(outname,'-djpeg','-r200');
%     close;
    end
%% plot coeff v.s. global mean T
elseif(choose_post == 2)
    for i_mode = 14

    load('pca_historical/Tmx_clim_historical_day.mat');
    % load tas
    load('global_stats/historical_tas_mean_1.mat');
    tas_historical = T_glb_mean;
    load('global_stats/historical_tas_mean_2.mat');
    tas_historical = [tas_historical; T_glb_mean];
    load('global_stats/historical_tas_mean_3.mat');
    tas_historical = [tas_historical; T_glb_mean];
    load('global_stats/historical_tas_mean_4.mat');
    tas_historical = [tas_historical; T_glb_mean];

    load('global_stats/ssp126_tas_mean_1.mat');
    tas_ssp126 = T_glb_mean;
    load('global_stats/ssp126_tas_mean_2.mat');
    tas_ssp126 = [tas_ssp126; T_glb_mean];

    load('global_stats/ssp585_tas_mean_1.mat');
    tas_ssp585 = T_glb_mean;
    load('global_stats/ssp585_tas_mean_2.mat');
    tas_ssp585 = [tas_ssp585; T_glb_mean];

    figure;
    Tm_historical = movmean(tas_historical - repmat(Tmx_clim,[165,1]),365);
    signal = movmean(coeff_historical(:,i_mode)/sqrt(lambda(i_mode)),365);
%     signal = coeff_historical(:,i_mode)/sqrt(lambda(i_mode));
    scatter(Tm_historical,signal,'k.','linewidth',1.2);
    hold on;
    Tm_ssp126 = movmean(tas_ssp126 - repmat(Tmx_clim,[86,1]),365);
    signal = movmean(coeff_ssp126(:,i_mode)/sqrt(lambda(i_mode)),365);
%     signal = coeff_ssp126(:,i_mode)/sqrt(lambda(i_mode));
    scatter(Tm_ssp126,signal,'g.','linewidth',1.2);
    Tm_ssp585 = movmean(tas_ssp585 - repmat(Tmx_clim,[86,1]),365);
    signal = movmean(coeff_ssp585(:,i_mode)/sqrt(lambda(i_mode)),365);
%     signal = coeff_ssp585(:,i_mode)/sqrt(lambda(i_mode));
    scatter(Tm_ssp585,signal,'r.','linewidth',1.2);

    % istart = (1979-1850)*365+1;
    % iend = (2019 - 1850 +1)*365;
    % signal_era5 = movmean(coeff_era5(:,i_mode)/sqrt(lambda(i_mode)),365);
    % plot(istart:iend,signal_era5,'Color',[0.5 0.5 0.5],'linewidth',1.2);

    xlim([min(Tm_historical),max(Tm_ssp585)]);
    box on;
    % ylim([-4 4]);
    set(gca,'TickLabelInterpreter', 'latex','fontsize',12);
    xlabel('$\overline{T}$','interpreter','latex','fontsize',20);
    ylabel('$\overline{a}_i / \sqrt{\lambda_i}$','interpreter','latex','fontsize',20);
%     set(gca,'XTick',ticks(1:20:end),...
%         'XTickLabel',ticklabel,'fontsize',12);
    set(gcf,'Units','inches',...
        'Position',[1 1 12 3], ...
        'PaperPosition',[0 0 12.5 3.5]);
%     legend('historical','ssp126','ssp585','location','northwest');
%     legend('boxoff');
%     outname = strcat('1011_pca_hist_coeff_Tm_',num2str(i_mode,'%3d'));
%     print(outname,'-djpeg','-r200');
%     close;
    end
end