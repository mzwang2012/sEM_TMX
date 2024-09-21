clear;
close all;

year_skip = 10;  % number of years sampled altogether



lat_skip = 1;
lon_skip = 1;
%% load grid and get the coarse grid
folder = 'download/';
variable = 'tasmax';

load('pca_historical/Tmx_clim_historical.mat','lon','lat','Tmx_mean');
n_lat = length(lat);
n_lon = length(lon);

n_latc = floor(n_lat/lat_skip);
n_lonc = floor(n_lon/lon_skip);

latc = lat(1:lat_skip:(n_latc-1)*lat_skip+1);
for i=1:lat_skip-1
    latc = latc + lat(1+i:lat_skip:(n_latc-1)*lat_skip+1+i);
end
latc = latc/lat_skip;

lonc = lon(1:lon_skip:(n_lonc-1)*lon_skip+1);
for i=1:lon_skip-1
    lonc = lonc + lon(1+i:lon_skip:(n_lonc-1)*lon_skip+1+i);
end
lonc = lonc/lon_skip;

theta = 90 - lat;
theta_rad = theta / 180 * pi;
phi_rad = lon / 180 *pi;
Mat_sin = reshape(sin(theta_rad),[1,n_lat]);


%% plot error v.s. t
load coastlines;
% n_ens = 25;
n_ens_list = [1,5,10,15,25,50];
% n_ens_cov_list = [10,15,20,25,50];
n_ens_cov_list = [10,15,25,50];

n_rlz = 25;
levels_Tmx = linspace(-2,12,20);
levels = linspace(-2,2,10);

% error_mean = zeros(25,1);


n = 1;
for k=1:6
n_ens = n_ens_list(k);
for year=2090 %1850:10:2090
    if(year < 2010)
        scenario = '';
    else
        scenario = 'ssp585_';
    end

    Tmx_std_true = zeros(n_lonc,n_latc);
    Tmx_mean_true = zeros(n_lonc,n_latc);
    for l=1:n_ens
        filename = strcat('stats_10year/ref/stats_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        Tmx_mean_true = Tmx_mean_true + TMX_avg;
        Tmx_std_true = Tmx_std_true + TMX_std.^2 + TMX_avg.^2;
    end
    Tmx_mean_true = Tmx_mean_true / n_ens;
    Tmx_std_true = Tmx_std_true / n_ens;
    Tmx_std_true = Tmx_std_true - Tmx_mean_true.^2;
    Tmx_std_true = sqrt(Tmx_std_true);

    q_true = zeros(n_lonc,n_latc);
    for l=1:n_ens

        filename = strcat('stats_10year/ref/qntl_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_true = q_true + qntl_975;
    end
    q_true = q_true / n_ens;


    Tmx_std_rom = zeros(size(Tmx_std_true));
    Tmx_std_rom_mrlz = zeros(n_lonc,n_latc,10);
    Tmx_mean_rom = zeros(size(Tmx_std_true));
    for l = 1:n_rlz
        filename = strcat(['stats_10year/rom_ens_',num2str(n_ens,'%d'), ...
            '/stats_rom_10year_'],scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        Tmx_mean_rom = Tmx_mean_rom + TMX_avg;
        Tmx_std_rom = Tmx_std_rom + TMX_std.^2 + TMX_avg.^2;
        Tmx_std_rom_mrlz(:,:,l) = TMX_std;
    end
    Tmx_mean_rom = Tmx_mean_rom / n_rlz;
    Tmx_std_rom = Tmx_std_rom / n_rlz;
    Tmx_std_rom = Tmx_std_rom - Tmx_mean_rom.^2;
    Tmx_std_rom = sqrt(Tmx_std_rom);

    q_rom_mrlz = zeros(n_lonc,n_latc,n_rlz);
    for l=1:n_rlz
        filename = strcat('stats_10year/rom_ens_',num2str(n_ens,'%d'), ...
            '/qntl_rom_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_rom_mrlz(:,:,l) = qntl_975;
    end
    q_rom = mean(q_rom_mrlz,3);

    % load COV ROM results
%     Tmx_std_ml = zeros(size(Tmx_std_true));
%     Tmx_std_ml_mrlz = zeros(n_lonc,n_latc,10);
%     Tmx_mean_ml = zeros(size(Tmx_std_true));
%     for l = 1:n_rlz
%         filename = strcat(['stats_10year/rom_cov/' ...
%             'stats_rom_cov_10year_'],scenario,num2str(year,'%d'),...
%             '_',num2str(l,'%03d'),'.mat');
%         load(filename);
%         Tmx_mean_ml = Tmx_mean_ml + TMX_avg;
%         Tmx_std_ml = Tmx_std_ml + TMX_std.^2 + TMX_avg.^2;
%         Tmx_std_ml_mrlz(:,:,l) = TMX_std;
%     end
%     Tmx_mean_ml = Tmx_mean_ml / n_rlz;
%     Tmx_std_ml = Tmx_std_ml / n_rlz;
%     Tmx_std_ml = Tmx_std_ml - Tmx_mean_ml.^2;
%     Tmx_std_ml = sqrt(Tmx_std_ml);

    error_mean_lat = trapz(phi_rad,(Tmx_mean_rom - Tmx_mean_true).^2,1);
    error_mean(n) = trapz(theta_rad,error_mean_lat.*Mat_sin,2) / (-4*pi);

    error_std_lat = trapz(phi_rad,(Tmx_std_rom - Tmx_std_true).^2,1);
    error_std(n) = trapz(theta_rad,error_std_lat.*Mat_sin,2) / (-4*pi);

    error_qntl_lat = trapz(phi_rad,(q_rom - q_true).^2,1);
    error_qntl(n) = trapz(theta_rad,error_qntl_lat.*Mat_sin,2) / (-4*pi);

    n = n + 1;
end
end

n = 1;
for k=1:length(n_ens_cov_list)
n_ens = n_ens_cov_list(k);
for year=2090 %1850:10:2090
    if(year < 2010)
        scenario = '';
    else
        scenario = 'ssp585_';
    end

    Tmx_std_true = zeros(n_lonc,n_latc);
    Tmx_mean_true = zeros(n_lonc,n_latc);
    for l=1:n_ens
        filename = strcat('stats_10year/ref/stats_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        Tmx_mean_true = Tmx_mean_true + TMX_avg;
        Tmx_std_true = Tmx_std_true + TMX_std.^2 + TMX_avg.^2;
    end
    Tmx_mean_true = Tmx_mean_true / n_ens;
    Tmx_std_true = Tmx_std_true / n_ens;
    Tmx_std_true = Tmx_std_true - Tmx_mean_true.^2;
    Tmx_std_true = sqrt(Tmx_std_true);

    q_true = zeros(n_lonc,n_latc);
    for l=1:n_ens

        filename = strcat('stats_10year/ref/qntl_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_true = q_true + qntl_975;
    end
    q_true = q_true / n_ens;


    Tmx_std_rom = zeros(size(Tmx_std_true));
    Tmx_std_rom_mrlz = zeros(n_lonc,n_latc,10);
    Tmx_mean_rom = zeros(size(Tmx_std_true));
    for l = 1:n_rlz
        filename = strcat(['stats_10year/rom_cov_ens_',num2str(n_ens,'%d'), ...
            '/stats_rom_cov_10year_'],scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        Tmx_mean_rom = Tmx_mean_rom + TMX_avg;
        Tmx_std_rom = Tmx_std_rom + TMX_std.^2 + TMX_avg.^2;
        Tmx_std_rom_mrlz(:,:,l) = TMX_std;
    end
    Tmx_mean_rom = Tmx_mean_rom / n_rlz;
    Tmx_std_rom = Tmx_std_rom / n_rlz;
    Tmx_std_rom = Tmx_std_rom - Tmx_mean_rom.^2;
    Tmx_std_rom = sqrt(Tmx_std_rom);

    q_rom_mrlz = zeros(n_lonc,n_latc,n_rlz);
    for l=1:n_rlz
        filename = strcat('stats_10year/rom_cov_ens_',num2str(n_ens,'%d'), ...
            '/qntl_rom_cov_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_rom_mrlz(:,:,l) = qntl_975;
    end
    q_rom = mean(q_rom_mrlz,3);

    % load COV ROM results
%     Tmx_std_ml = zeros(size(Tmx_std_true));
%     Tmx_std_ml_mrlz = zeros(n_lonc,n_latc,10);
%     Tmx_mean_ml = zeros(size(Tmx_std_true));
%     for l = 1:n_rlz
%         filename = strcat(['stats_10year/rom_cov/' ...
%             'stats_rom_cov_10year_'],scenario,num2str(year,'%d'),...
%             '_',num2str(l,'%03d'),'.mat');
%         load(filename);
%         Tmx_mean_ml = Tmx_mean_ml + TMX_avg;
%         Tmx_std_ml = Tmx_std_ml + TMX_std.^2 + TMX_avg.^2;
%         Tmx_std_ml_mrlz(:,:,l) = TMX_std;
%     end
%     Tmx_mean_ml = Tmx_mean_ml / n_rlz;
%     Tmx_std_ml = Tmx_std_ml / n_rlz;
%     Tmx_std_ml = Tmx_std_ml - Tmx_mean_ml.^2;
%     Tmx_std_ml = sqrt(Tmx_std_ml);

    error_cov_mean_lat = trapz(phi_rad,(Tmx_mean_rom - Tmx_mean_true).^2,1);
    error_cov_mean(n) = trapz(theta_rad,error_cov_mean_lat.*Mat_sin,2) / (-4*pi);

    error_cov_std_lat = trapz(phi_rad,(Tmx_std_rom - Tmx_std_true).^2,1);
    error_cov_std(n) = trapz(theta_rad,error_cov_std_lat.*Mat_sin,2) / (-4*pi);

    error_cov_qntl_lat = trapz(phi_rad,(q_rom - q_true).^2,1);
    error_cov_qntl(n) = trapz(theta_rad,error_cov_qntl_lat.*Mat_sin,2) / (-4*pi);

    n = n + 1;
end
end

error_mean = sqrt(error_mean);
error_std = sqrt(error_std);
error_qntl = sqrt(error_qntl);
error_cov_mean = sqrt(error_cov_mean);
error_cov_std = sqrt(error_cov_std);
error_cov_qntl = sqrt(error_cov_qntl);
%% plot
figure;
plot(n_ens_list,error_mean,'bo--','linewidth',1.2);
hold on;
plot(n_ens_list,error_std,'go--','linewidth',1.2);
plot(n_ens_list,error_qntl,'ko--','linewidth',1.2);

plot(n_ens_cov_list,error_cov_mean,'bo-','linewidth',1.2);
hold on;
plot(n_ens_cov_list,error_cov_std,'go-','linewidth',1.2);
plot(n_ens_cov_list,error_cov_qntl,'ko-','linewidth',1.2);
xlim([0 50]);
% set(gca,'XScale','Log');
set(gca,'fontsize',12,'TickLabelInterpreter','latex');
xlabel('$N_{\omega}$','interpreter','latex','fontsize',15);
ylabel('RMSE','interpreter','latex','fontsize',15);
set(gcf,'Units','inches','Position',[1,1,4.2,3]);

figure;
plot(n_ens_list,error_mean/error_mean(1),'bo--','linewidth',1.2);
hold on;
plot(n_ens_list,error_std/error_std(1),'go--','linewidth',1.2);
plot(n_ens_list,error_qntl/error_qntl(1),'ko--','linewidth',1.2);

plot(n_ens_cov_list,error_cov_mean/error_cov_mean(1),'bo-','linewidth',1.2);
hold on;
plot(n_ens_cov_list,error_cov_std/error_cov_std(1),'go-','linewidth',1.2);
plot(n_ens_cov_list,error_cov_qntl/error_cov_qntl(1),'ko-','linewidth',1.2);
xlim([0 50]);
% set(gca,'XScale','Log');
set(gca,'fontsize',12,'TickLabelInterpreter','latex');
xlabel('$N_{\omega}$','interpreter','latex','fontsize',15);
ylabel('$\textrm{RMSE}/\textrm{RMSE}_0$','interpreter','latex','fontsize',15);
set(gcf,'Units','inches','Position',[1,1,4.2,3]);