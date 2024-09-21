clear;
close all;

choose_data = 0; % 1 for loading data from mounted folder; 
choose_post = 1; % 1 for TMX averaged over 10 years
                 % 2 for std of TMX
                 % 6 for quantiles, 10 years
year_skip = 10;  % number of years sampled altogether

lat_skip = 2;
lon_skip = 2;
%% load grid and get the coarse grid
folder = 'download/';
variable = 'tasmax';

filename = 'pca_ERA5/Tmx_365_global_1979_1998_modes.mat';
load(filename,'lat','lon');

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

%% post processing
startYear = 2010;
endYear = 2099;

n_rlz = 10;
n_10years = (endYear - startYear + 1)/10;
true_mean = zeros(n_10years,4);
true_std = zeros(n_10years,4);
true_qntl = zeros(n_10years,4);
error_mean = zeros(n_10years,4);
error_std = zeros(n_10years,4);
error_qntl = zeros(n_10years,4);
folder_stats = ["stats_winter/","stats_spring/","stats_10year/","stats_autumn/"];

theta = 90 - latc;
theta_rad = theta / 180 * pi;
phi_rad = lonc / 180 *pi;
Mat_sin = reshape(sin(theta_rad),[1,n_latc]);

for s=1:4
    n = 1;
    for year = startYear:year_skip:endYear

        if(year < 2009)
            scenario = '';
        else
            scenario = 'ssp585_';
        end
        filename = strcat(folder_stats(s),'stats_10year_',scenario,num2str(year,'%d'),'.mat');
        load(filename);
        Tmx_mean_true = TMX_avg;
        Tmx_std_true = TMX_std;

        Tmx_std_rom = zeros(size(Tmx_std_true));
        Tmx_std_rom_mrlz = zeros(n_lonc,n_latc,10);
        Tmx_mean_rom = zeros(size(Tmx_std_true));
        for l = 1:n_rlz
            filename = strcat(folder_stats(s),'stats_rom_10year_',scenario,num2str(year,'%d'),...
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

        % load ML model results
        %     Tmx_std_ml = zeros(size(Tmx_std_true));
        %     Tmx_std_ml_mrlz = zeros(n_lonc,n_latc,10);
        %     Tmx_mean_ml = zeros(size(Tmx_std_true));
        %     for l = 1:n_rlz
        %         filename = strcat(folder_stats,'stats_rom_10year_',scenario,num2str(year,'%d'),...
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

        filename = strcat(folder_stats(s),'qntl_10year_',scenario,num2str(year,'%d'),'.mat');
        load(filename);
        q_true = qntl_975;

        q_rom_mrlz = zeros(n_lonc,n_latc,n_rlz);
        for l=1:n_rlz
            filename = strcat(folder_stats(s),'qntl_rom_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            q_rom_mrlz(:,:,l) = qntl_975;
        end
        q_rom = mean(q_rom_mrlz,3);

        %     q_ml_mrlz = zeros(n_lonc,n_latc,n_rlz);
        %     for l=1:n_rlz
        %         filename = strcat(folder_stats,'qntl_rom_10year_',scenario,num2str(year,'%d'),...
        %             '_',num2str(l,'%03d'),'.mat');
        %         load(filename);
        %         q_ml_mrlz(:,:,l) = qntl_975;
        %     end
        %     q_ml = mean(q_ml_mrlz,3);

        true_mean_lat = trapz(phi_rad,Tmx_mean_true,1);
        true_mean(n,s) = trapz(theta_rad,true_mean_lat.*Mat_sin,2) / (-4*pi);

        true_std_lat = trapz(phi_rad,Tmx_std_true.^2,1);
        true_std(n,s) = trapz(theta_rad,true_std_lat.*Mat_sin,2) / (-4*pi);

        true_qntl_lat = trapz(phi_rad,q_true,1);
        true_qntl(n,s) = trapz(theta_rad,true_qntl_lat.*Mat_sin,2) / (-4*pi);

        error_mean_lat = trapz(phi_rad,(Tmx_mean_rom - Tmx_mean_true).^2,1);
        error_mean(n,s) = trapz(theta_rad,error_mean_lat.*Mat_sin,2) / (-4*pi);

        error_std_lat = trapz(phi_rad,(Tmx_std_rom - Tmx_std_true).^2,1);
        error_std(n,s) = trapz(theta_rad,error_std_lat.*Mat_sin,2) / (-4*pi);

        error_qntl_lat = trapz(phi_rad,(q_rom - q_true).^2,1);
        error_qntl(n,s) = trapz(theta_rad,error_qntl_lat.*Mat_sin,2) / (-4*pi);

        n = n + 1;
    end
end

true_std = sqrt(true_std);
error_mean = sqrt(error_mean);
error_std = sqrt(error_std);
error_qntl = sqrt(error_qntl);

%%
error_ssp126_mean = error_mean;
error_ssp126_std = error_std;
error_ssp126_qntl = error_qntl;
load('rmse_ssp585.mat');
time_ssp585 = 1850:10:2090;
time_ssp126 = 2010:10:2090;

figure;
subplot(3,1,1);
plot(time_ssp585,error_ssp585_mean(:,1),'b','linewidth',1.5);
hold on;
plot(time_ssp585,error_ssp585_mean(:,2),'g','linewidth',1.5);
plot(time_ssp585,error_ssp585_mean(:,3),'r','linewidth',1.5);
plot(time_ssp585,error_ssp585_mean(:,4),'y','linewidth',1.5);

plot(time_ssp126,error_ssp126_mean(:,1),'b--','linewidth',1.5);
hold on;
plot(time_ssp126,error_ssp126_mean(:,2),'g--','linewidth',1.5);
plot(time_ssp126,error_ssp126_mean(:,3),'r--','linewidth',1.5);
plot(time_ssp126,error_ssp126_mean(:,4),'y--','linewidth',1.5);

xlim([1850 2090]);
ylim([0.2 1]);
xlabel('year','interpreter','latex','fontsize',15);
ylabel('RMSE of mean','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter', 'latex','fontsize',15);

h = legend('Dec-Feb','Mar-May','Jun-Aug','Sep-Nov',...
    'orientation','horizontal','location','northeast');
set(h,'box','off');
set(h,'interpreter','latex');
% annotation('line',[0.7,0.75],[0.8,0.8],'Color','k','linestyle','-');
% text(0.76,0.79,'historical, SSP5-8.5','fontsize',10);

subplot(3,1,2);
plot(time_ssp585,error_ssp585_std(:,1),'b','linewidth',1.5);
hold on;
plot(time_ssp585,error_ssp585_std(:,2),'g','linewidth',1.5);
plot(time_ssp585,error_ssp585_std(:,3),'r','linewidth',1.5);
plot(time_ssp585,error_ssp585_std(:,4),'y','linewidth',1.5);

plot(time_ssp126,error_ssp126_std(:,1),'b--','linewidth',1.5);
hold on;
plot(time_ssp126,error_ssp126_std(:,2),'g--','linewidth',1.5);
plot(time_ssp126,error_ssp126_std(:,3),'r--','linewidth',1.5);
plot(time_ssp126,error_ssp126_std(:,4),'y--','linewidth',1.5);

xlim([1850 2090]);
% ylim([0.2 0.8]);
xlabel('year','interpreter','latex','fontsize',15);
ylabel('RMSE of std','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter', 'latex','fontsize',15);


subplot(3,1,3);
plot(time_ssp585,error_ssp585_qntl(:,1),'b','linewidth',1.5);
hold on;
plot(time_ssp585,error_ssp585_qntl(:,2),'g','linewidth',1.5);
plot(time_ssp585,error_ssp585_qntl(:,3),'r','linewidth',1.5);
plot(time_ssp585,error_ssp585_qntl(:,4),'y','linewidth',1.5);

plot(time_ssp126,error_ssp126_qntl(:,1),'b--','linewidth',1.5);
hold on;
plot(time_ssp126,error_ssp126_qntl(:,2),'g--','linewidth',1.5);
plot(time_ssp126,error_ssp126_qntl(:,3),'r--','linewidth',1.5);
plot(time_ssp126,error_ssp126_qntl(:,4),'y--','linewidth',1.5);

xlim([1850 2090]);
ylim([0.5 3]);
xlabel('year','interpreter','latex','fontsize',15);
ylabel('RMSE of qntl','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter', 'latex','fontsize',15);
set(gcf,'Units','inches',...
    'Position',[1 1 12 7]);
