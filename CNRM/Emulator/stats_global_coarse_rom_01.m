clear;
close all;

%% path to CMIP6 data
% folder = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/historical/';
% folder = 'download/';
folder = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/';
scenario = 'ssp585';
variable = 'tasmax';
flag_season = 4; % 1 for winter, 2 for spring, 3 for summer, 4 for autumn

currentFile = mfilename;
fileIdx = str2num(currentFile(end-1:end));
disp(fileIdx);
n_seed = fileIdx;
rng(n_seed);

skipYear = 10;  % number of years sampled altogether
n_eig = 50*40;

lat_skip = 2;
lon_skip = 2;

% load coordinates
filename = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/historical/tasmax/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_18500101-18991231.nc';
%filename = 'download/historical/tasmax/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_18500101-18991231.nc';
lat = ncread(filename,'lat');
n_lat = length(lat);
lon = ncread(filename,'lon');
n_lon = length(lon);

%% load global stats
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

%% load climatological mean and PCA modes
inname = strcat('pca_historical/modes_Tmx_365skip5_global_historical.mat');
load(inname,'lambda','w_hat');
mode = w_hat;
clear w_hat;
load(inname,'lambda','w_hat_2');
mode(:,:,501:1000) = w_hat_2;
clear w_hat_2;

if(n_eig > 1000)
    inname = strcat('pca_historical/modes2_Tmx_365skip5_global_historical.mat');
    load(inname,'lambda','w_hat_3');
    mode(:,:,1001:1500) = w_hat_3;
    clear w_hat_3;
    load(inname,'lambda','w_hat_4');
    mode(:,:,1501:2000) = w_hat_4;
    clear w_hat_4;
end

[n_lon,n_lat,~] = size(mode);

load('pca_historical/Tmx_clim_historical.mat','Tmx_mean');

startYear = 1850;
endYear = 2009;
n_years = endYear - startYear + 1;
n_days =  n_years * 365;

% load coefficients
coeff_historical = [];
for k = 1:(n_eig/50)
    filename = strcat('pca_historical/Tmx_365_global_historical_coeff_',num2str(k,'%02d'),'.mat');
    load(filename,'coeff');
    coeff_historical = [coeff_historical,coeff];
end

% use ssp 585 for training
coeff_ssp = [];
for k = 1:(n_eig/50)
    filename = strcat('pca_historical/Tmx_365_global_ssp585_coeff_',num2str(k,'%02d'),'.mat');
    load(filename,'coeff');
    coeff_ssp = [coeff_ssp,coeff];
end

%% ROM
% modeling
data_est = [coeff_historical; coeff_ssp];
T_glb = [tas_historical;tas_ssp585];
T_glb = T_glb - repmat(Tmx_clim,[251,1]);
T_glb_ssp126 = [tas_historical;tas_ssp126];
T_glb_ssp126 = T_glb_ssp126 - repmat(Tmx_clim,[251,1]);
t_est = 1:(365*251);
%mu_pred_all = [];
%var_pred_all = [];
%mu_valid_all = [];
%var_valid_all = [];
%for k=1:40
%    filename = strcat('ML_output/mu_var_valid_',num2str(k,'%02d'));
%    load(filename);
%    mu_pred_all = [mu_pred_all,mu_pred];
%    var_pred_all = [var_pred_all,var_pred];
%    mu_valid_all = [mu_valid_all,mu_valid];
%    var_valid_all = [var_valid_all,var_valid];
%end
%% coeff_rom = rom_glbT_ML_mean_var_pred(data_est,t_est,t_est,mu_pred_all,var_pred_all,mu_valid_all,var_valid_all);
%if(scenario == 'ssp585')
%    coeff_rom = rom_glbT_ML_mean_var_pred(data_est,t_est,t_est,mu_pred_all,var_pred_all,mu_pred_all,var_pred_all);
%elseif(scenario == 'ssp126')
%    coeff_rom = rom_glbT_ML_mean_var_pred(data_est,t_est,t_est,mu_pred_all,var_pred_all,mu_valid_all,var_valid_all);
%end

if(scenario == 'ssp585')
    coeff_rom = rom_glbT_lin_mean_var_pred(data_est,T_glb,T_glb,t_est,t_est);
elseif(scenario == 'ssp126')
    coeff_rom = rom_glbT_lin_mean_var_pred(data_est,T_glb,T_glb_ssp126,t_est,t_est);
end

% Going to average over a local neighborhood
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

switch flag_season
    case 1
        folder_output = 'stats_winter/';
        day_start = 1;
        day_end = 59;
        ndays_season = 90;
        nlength = 90*lat_skip*lon_skip;
    case 2
        folder_output = 'stats_spring/';
        day_start = 60;
        day_end = 151;
        nlength = 92*lat_skip*lon_skip;
        ndays_season = 92;
    case 3
        folder_output = 'stats_summer/';
        day_start = 152;
        day_end = 243;
        nlength = 92*lat_skip*lon_skip;
        ndays_season = 92;
    case 4
        folder_output = 'stats_autumn/';
        day_start = 244;
        day_end = 334;
        nlength = 91*lat_skip*lon_skip;
        ndays_season = 91;
    otherwise
        disp('flag_season must be 1, 2, 3, 4');
end

%% post processing
for statsYear = startYear:skipYear:endYear
    currentDate = datetime(statsYear,1,1);

    % collect samples from the local neighborhood
    Tmx_c = zeros(n_lonc,n_latc,ndays_season*lat_skip*lon_skip*skipYear);
    nn = 0;
    for readYear = statsYear:statsYear+skipYear-1    
        disp(readYear);
        % PCA reconstruction
        Tmx_est = zeros(n_lon,n_lat,365);
        %if(readYear < 2015)
        k = readYear - 1850 + 1;
        tic;
        for j=day_start:day_end
        for i_mode = 1:n_eig
            indx = 365*(k-1)+j;
            Tmx_est(:,:,j) = Tmx_est(:,:,j) + coeff_rom(indx,i_mode) * mode(:,:,i_mode);
        end
        end
        if(flag_season == 1)
            for j=335:365
            for i_mode = 1:n_eig
                indx = 365*(k-1)+j;
                Tmx_est(:,:,j) = Tmx_est(:,:,j) + coeff_rom(indx,i_mode) * mode(:,:,i_mode);
            end
            end
        end
        toc;
        %else
        %    k = readYear - 2015 + 1;
        %    for j=day_start:day_end
        %    for i_mode = 1:n_eig
        %        indx = 365*(k-1)+j;
        %        Tmx_est(:,:,j) = Tmx_est(:,:,j) + coeff_ssp(indx,i_mode) * mode(:,:,i_mode);
        %    end
        %    end
        %end
        % add climatological mean
        Tmx_est(:,:,day_start:day_end) = Tmx_est(:,:,day_start:day_end) + ...
                                     Tmx_mean(:,:,day_start:day_end);
        if(flag_season == 1)
            Tmx_est(:,:,335:365) = Tmx_est(:,:,335:365) + ...
                                         Tmx_mean(:,:,335:365);
	end
        % remove zero elements
        if(flag_season == 1)
	    Tmx = zeros(n_lon,n_lat,90);
            Tmx(:,:,1:59) = Tmx_est(:,:,1:59);
            Tmx(:,:,60:90) = Tmx_est(:,:,335:365);

        else
            Tmx = Tmx_est(:,:,day_start:day_end);
        end
 
        clear Tmx_est;
        tic;
        for j=1:n_latc
            for i=1:n_lonc
                istart = (i-1)*lon_skip+1;
                iend = i*lon_skip;
                jstart = (j-1)*lat_skip+1;
                jend = j*lat_skip;
                tmp = Tmx(istart:iend,jstart:jend,:);
                Tmx_c(i,j,nlength*nn+1:nlength*(nn+1)) = reshape(tmp,[1,1,nlength]);
            end
        end
        toc;
        
        nn = nn + 1;
     end

     % mean and standard deviation
     TMX_avg = mean(Tmx_c,3);
     TMX_std = std(Tmx_c,0,3);
     
     if(statsYear + skipYear - 1 < 2015)
         outname = strcat(folder_output,'stats_rom_10year_',num2str(statsYear,'%d'),'_',num2str(n_seed,'%03d'),'.mat');
     else
         outname = strcat(folder_output,'stats_rom_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(n_seed,'%03d'),'.mat');
     end
     save(outname,'TMX_avg','TMX_std','lonc','latc','skipYear');

     %% compute quantiles
     qntl_75 = zeros(n_lonc,n_latc);
     qntl_90 = zeros(n_lonc,n_latc);
     qntl_95 = zeros(n_lonc,n_latc);
     qntl_975 = zeros(n_lonc,n_latc);
     for j=1:n_latc
         for i=1:n_lonc
             qntl_75(i,j) = quantile(Tmx_c(i,j,:),0.75);
             qntl_90(i,j) = quantile(Tmx_c(i,j,:),0.9);
             qntl_95(i,j) = quantile(Tmx_c(i,j,:),0.95);
             qntl_975(i,j) = quantile(Tmx_c(i,j,:),0.975);
         end
     end
     if(statsYear + skipYear - 1 < 2015)
         outname = strcat(folder_output,'qntl_rom_10year_',num2str(statsYear,'%d'),'_',num2str(n_seed,'%03d'),'.mat');
     else
         outname = strcat(folder_output,'qntl_rom_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(n_seed,'%03d'),'.mat');
     end
     save(outname,'qntl_75','qntl_90','qntl_95','qntl_975');

    %% extract local PDF
     n_samp = ndays_season*lat_skip*lon_skip*skipYear;
     Tmx_local = zeros(n_samp,19);
     city_coords = readtable('stats_summer/city_coords.csv');
     for i=1:19
         if(city_coords.lon(i) < 0)
             lon_tmp = city_coords.lon(i) + 360;
         else
             lon_tmp = city_coords.lon(i);
         end
         if(lon_tmp > 359.5)
             lon_tmp = 0;
         end
         idx_lon = find(abs(lonc - lon_tmp) < 0.250*lon_skip);
         idx_lat = find(abs(latc - city_coords.lat(i)) < 0.250*lat_skip);
         disp([lonc(idx_lon),latc(idx_lat)]);
         tmp = reshape(Tmx_c(idx_lon,idx_lat,:),[n_samp,1]);
         Tmx_local(:,i) = tmp; 
     end

     if(statsYear + skipYear - 1 < 2015)
         outname = strcat(folder_output,'localPDF_rom_10year_',num2str(statsYear,'%d'),'_',num2str(n_seed,'%03d'),'.mat');;
     else
         outname = strcat(folder_output,'localPDF_rom_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(n_seed,'%03d'),'.mat');;
     end
     save(outname,'Tmx_local','city_coords');
end
