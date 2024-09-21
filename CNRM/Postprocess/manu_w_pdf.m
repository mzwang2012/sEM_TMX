clear;
close all;

choose_post = 1; % 0 for loading coefficients, run ROM to get w
                 % 1 for loading w

n_eig = 50*40;
n_rlz = 1;
rng('default');

% rng(1);
%% load coefficients
if(choose_post == 0)
    load('pca_historical/modes_Tmx_365skip5_global_historical.mat','lambda');
    coeff_historical = [];
    % for k=38:38+(n_eig/50)
    for k=1:(n_eig/50)
        filename = strcat('pca_historical/Tmx_365_global_historical_coeff_',num2str(k,'%02d'),'.mat');
        load(filename);
        coeff_historical = [coeff_historical,coeff];
    end
    % coeff_historical = coeff;

    coeff_ssp585 = [];
    % for k=38:38+(n_eig/50)
    for k=1:(n_eig/50)
        filename = strcat('pca_historical/Tmx_365_global_ssp585_coeff_',num2str(k,'%02d'),'.mat');
        load(filename);
        coeff_ssp585 = [coeff_ssp585,coeff];
    end

    % coeff_ssp126 = [];
    % % for k=38:38+(n_eig/50)
    % for k=1:(n_eig/50)
    %     filename = strcat('pca_historical/Tmx_365_global_ssp126_coeff_',num2str(k,'%02d'),'.mat');
    %     load(filename);
    %     coeff_ssp126 = [coeff_ssp126,coeff];
    % end

    disp('finished loading coeffs');

    % load('../ERA5/PCA/Tmx_365_global_1979_2019_coeff.mat');
    % coeff_era5 = coeff;
    % clear coeff;

    n_days_historical = (2015 - 1850)*365;
    time_historical = 1:n_days_historical;
    n_days_ssp = (2100-2015+1)*365;
    time_ssp = (n_days_historical+1):(n_days_historical+n_days_ssp);

    lambda = lambda / (n_days_historical/5);

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

    data_est = [coeff_historical(:,1:n_eig); coeff_ssp585(:,1:n_eig)];
    T_glb = [tas_historical;tas_ssp585];
    T_glb = T_glb - repmat(Tmx_clim,[251,1]);
    t_est = 1:(n_days_historical+n_days_ssp);
    [true_w,rom_w] = rom_glbT_lin_mean_var_w(data_est,T_glb,T_glb,t_est,t_est);

elseif(choose_post == 1)
    %%
    load('true_rom_w.mat','true_w');
    true_w_summer = true_w{3};
    clear true_w;
    load('rom_w_summer.mat');

    %%
    xgrid = -5:0.01:5;

    figure;

    % first panel
    i_mode = 1;
    [f,xi] = ksdensity(true_w_summer(:,i_mode));
    data_rom = reshape(rom_w_summer_mrlz(:,i_mode,:),[23092*10,1]);
    f_rom = ksdensity(data_rom,xi);

    [mu,sigma] = normfit(true_w_summer(:,i_mode));
    f_normal = 1/(sigma*sqrt(2*pi)) * exp(-(xgrid - mu).^2/(2*sigma^2));
    
    subplot(1,3,1);
    plot(xgrid,f_normal,'Color',[0.5 0.5 0.5,0.7],'linewidth',5);
    hold on;
    plot(xi,f_rom,'b','linewidth',1.5);
    hold on;
    plot(xi,f,'r--','linewidth',1.5);
    axis([-5 5 0 0.42]);
    set(gca,'fontsize',13,'TickLabelInterpreter','latex');
    xlabel('$\eta_{s,i}$','fontsize',15,'interpreter','latex');
    ylabel('PDF','fontsize',13);
    text(-0.22,0.95,'($a$)','interpreter','latex','Units','normalized','FontSize',14);
    h = legend('Gaussian','Emulation','Data','location','northwest');
    set(h,'box','off','fontsize',9);
    

    % second panel
    i_mode = 2;
    [f,xi] = ksdensity(true_w_summer(:,i_mode));
    data_rom = reshape(rom_w_summer_mrlz(:,i_mode,:),[23092*10,1]);
    f_rom = ksdensity(data_rom,xi);

    [mu,sigma] = normfit(true_w_summer(:,i_mode));
    f_normal = 1/(sigma*sqrt(2*pi)) * exp(-(xgrid - mu).^2/(2*sigma^2));
    
    subplot(1,3,2);
    plot(xgrid,f_normal,'Color',[0.5 0.5 0.5,0.7],'linewidth',5);
    hold on;
    plot(xi,f_rom,'b','linewidth',1.5);
    hold on;
    plot(xi,f,'r--','linewidth',1.5);
    axis([-5 5 0 0.42]);
    set(gca,'fontsize',13,'TickLabelInterpreter','latex');
    xlabel('$\eta_{s,i}$','fontsize',15,'interpreter','latex');
    ylabel('PDF','fontsize',13);
    text(-0.22,0.95,'($b$)','interpreter','latex','Units','normalized','FontSize',14)
    
    % third panel
    i_mode = 500;
    [f,xi] = ksdensity(true_w_summer(:,i_mode));
    data_rom = reshape(rom_w_summer_mrlz(:,i_mode,:),[23092*10,1]);
    f_rom = ksdensity(data_rom,xi);

    [mu,sigma] = normfit(true_w_summer(:,i_mode));
    f_normal = 1/(sigma*sqrt(2*pi)) * exp(-(xgrid - mu).^2/(2*sigma^2));
    
    subplot(1,3,3);
    plot(xgrid,f_normal,'Color',[0.5 0.5 0.5,0.7],'linewidth',5);
    hold on;
    plot(xi,f_rom,'b','linewidth',1.5);
    hold on;
    plot(xi,f,'r--','linewidth',1.5);
    axis([-5 5 0 0.42]);
    set(gca,'fontsize',13,'TickLabelInterpreter','latex');
    xlabel('$\eta_{s,i}$','fontsize',15,'interpreter','latex');
    ylabel('PDF','fontsize',13);
    text(-0.22,0.95,'($c$)','interpreter','latex','Units','normalized','FontSize',14)
    

    set(gcf,'Units','inches',...
        'Position',[1 1 13 3]);

end