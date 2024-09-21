clear;
close all;

%currentFile = mfilename;
%fileIdx = str2num(currentFile(end-1:end));
%disp(fileIdx);
%n_seed = fileIdx;
rng(1);

%% path to CMIP6 data
folder = '/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/';
scenario = 'ssp585';
variable = 'tasmax';

n_eig = 1000;

mem_start = 1; 
mem_end = 15;
n_ens = (mem_end - mem_start)+1;
n_rlz = 50;

lat_skip = 1;
lon_skip = 1;
%% load climatological mean and PCA modes
load('pca_historical/Tmx_clim_historical.mat','lon','lat','Tmx_mean');
n_lat = length(lat);
n_lon = length(lon);

% load('pca_ERA5/Tmx_clim_1979_1998.mat','Tmx_mean');

startYear = 1850; %1850;
endYear = 2100; %2014;
skipYear = 10;  % number of years sampled altogether
n_years = endYear - startYear + 1;
n_days =  n_years * 365;

n_latc = floor(n_lat/lat_skip);
n_lonc = floor(n_lon/lon_skip);

% Going to average over a local neighborhood
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

%% load global stats
load('global_stats/historical_tas_mean.mat','T_glb_mean');
tas_historical = T_glb_mean;
load('global_stats/ssp585_tas_mean.mat','T_glb_mean');
tas_ssp585 = T_glb_mean;
load('global_stats/ssp119_tas_mean.mat','T_glb_mean');
tas_ssp119 = T_glb_mean;


%% load PCA coefficients
disp('start loading pca coefficients');
l1 = length(tas_historical) + length(tas_ssp585);
coeff_ssp = zeros(l1,n_eig,n_ens);
for k = 1:(n_eig/50)
    disp(k);
    filename = strcat('pca_historical/Tmx_365_global_ssp585_coeff_',num2str(k,'%02d'),'.mat');
    %% DANGER ZONE!!!!!!!!!! Always use ssp585 for training !!!!!!!!!! %%
    load(filename,'coeff');
    istart = (k-1)*50+1;
    iend = k*50;
    coeff_ssp(:,istart:iend,:) = coeff(:,:,mem_start:mem_end);
end
disp('finish loading pca coefficients');

%% ROM
% modeling
T_glb = [tas_historical;tas_ssp585];
T_glb = T_glb - 273.15;
T_glb_ssp119 = [tas_historical;tas_ssp119];
T_glb_ssp119 = T_glb_ssp119 - 273.15;
%T_glb = T_glb - repmat(Tmx_clim,[251,1]);
%T_glb_ssp126 = [tas_historical;tas_ssp126];
%T_glb_ssp126 = T_glb_ssp126 - repmat(Tmx_clim,[251,1]);
t_est = 1:(365*251);
coeff_rom_full = zeros(l1,n_eig,n_rlz);
for i=1:5
    j_start = (i-1)*10+1;
    j_end = i*10;
    coeff_rom_full(:,:,j_start:j_end) = rom_LE_glbT_lin_mean_cov_pred(coeff_ssp(:,:,mem_start:mem_end),T_glb,T_glb,t_est,t_est,j_end-j_start+1);
end
%coeff_rom_full = rom_LE_glbT_lin_mean_cov_pred(coeff_ssp(:,:,mem_start:mem_end),T_glb,T_glb,t_est,t_est,n_rlz);
clear coeff_ssp;


%% output
daystart = (startYear - 1850)*365 + 1;
dayend   = (endYear - 1850 + 1)*365;
disp([daystart,dayend]);
for k=1:(n_eig/50)
    istart = (k-1)*50+1;
    iend = k*50;
    coeff_rom = coeff_rom_full(daystart:dayend,istart:iend,:);
    filename = strcat('rom_output/coeff_rom_cov_',scenario,'_ens_',num2str(n_ens,'%02d'),'_',num2str(k,'%02d'),'.mat');
    save(filename,'coeff_rom','startYear','endYear','mem_start','mem_end');
end
