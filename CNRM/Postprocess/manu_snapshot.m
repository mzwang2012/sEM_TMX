clear;
close all;

folder = 'download/';
variable = 'tasmax';

filename = 'pca_ERA5/Tmx_365_global_1979_1998_modes.mat';
load(filename,'lat','lon');
n_lat = length(lat);
n_lon = length(lon);

Tmx = ncread('mnt/ssp585/tasmax/tasmax_day_CNRM-CM6-1-HR_ssp585_r1i1p1f2_gr_20650101-21001231.nc',...
    'tasmax',[1,1,12934],[n_lon,n_lat,1]);

load('stats_noice/stats_historical_1850_1899.mat','Tmx_mean_ref','Tmx_std_ref','q_ref');


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

%%
load('pca_noice/ice_mask.mat');
no_ice = ~(is_ice');
no_ice = [no_ice(n_lon/2+1:n_lon,:); no_ice(1:n_lon/2,:)];
% idx = ~is_ice;
% Tmx(~no_ice) = NaN;

Tmx = Tmx(1:2:end,1:2:end) - Tmx_mean_ref;

lon_plot = [lonc(n_lonc/2+1:n_lonc) - 360;lonc(1:n_lonc/2)];
mu_plot = [Tmx(n_lonc/2+1:n_lonc,:); Tmx(1:n_lonc/2,:)];


load coastlines;
figure;
% set(gcf,'Units','inches',...
%     'Position',[1 1 20 6],'PaperPosition',[0 0 21 7]);

% h1 = subplot('Position',[0.05 0.05 0.27 0.9]);
ax = worldmap('world');
contourfm(latc,lon_plot,mu_plot',20,'linestyle','none');
hold on;
geoshow(coastlat,coastlon,'color','k');
cb = colorbar;
colormap(flipud(hot));
caxis([-2 12]);
set(gca,'fontsize',12);
% setm(ax,'FFaceColor', [0.7 0.7 0.7]);
mlabel('off');
set(gca,'TickLabelInterpreter', 'latex');