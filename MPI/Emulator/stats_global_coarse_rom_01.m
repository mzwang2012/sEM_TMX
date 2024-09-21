clear;
close all;

%currentFile = mfilename;
%fileIdx = str2num(currentFile(end-2:end));
%disp(fileIdx);
%
%mem_id = fileIdx;

currentFile = mfilename;
fileIdx = str2num(currentFile(end-1:end));
disp(fileIdx);
n_seed = fileIdx;
rng(n_seed);

%% path to CMIP6 data
folder = '/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/';
folder_coeff = '/net/fs06/d3/mzwang/ROM_MPI_LR_tasmax/rom_output/';
file_coeff = 'coeff_rom_';
scenario = 'ssp126';
variable = 'tasmax';

n_eig = 1000;

mem_start = 1; 
mem_end = 50;
n_ens = (mem_end - mem_start)+1;
folder_output = strcat('stats_10year/rom_ens_',num2str(n_ens,'%d'),'/');
file_output = 'rom';

lat_skip = 1;
lon_skip = 1;
%% load climatological mean and PCA modes
load('pca_historical/Tmx_clim_historical.mat','lon','lat','Tmx_mean');
n_lat = length(lat);
n_lon = length(lon);

% load('pca_ERA5/Tmx_clim_1979_1998.mat','Tmx_mean');

startYear = 2090; %1850;
endYear = 2099; %2014;
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
%disp('start loading pca coefficients');
l1 = length(tas_historical) + length(tas_ssp585);
%coeff_ssp = zeros(l1,n_eig,n_ens);
%for k = 1:(n_eig/50)
%    disp(k);
%    filename = strcat('pca_historical/Tmx_365_global_ssp585_coeff_',num2str(k,'%02d'),'.mat');
%    %% DANGER ZONE!!!!!!!!!! Always use ssp585 for training !!!!!!!!!! %%
%    load(filename,'coeff');
%    istart = (k-1)*50+1;
%    iend = k*50;
%    coeff_ssp(:,istart:iend,:) = coeff(:,:,mem_start:mem_end);
%end
%disp('finish loading pca coefficients');
%
%%% ROM
%% modeling
%T_glb = [tas_historical;tas_ssp585];
%T_glb = T_glb - 273.15;
%T_glb_ssp119 = [tas_historical;tas_ssp119];
%T_glb_ssp119 = T_glb_ssp119 - 273.15;
%%T_glb = T_glb - repmat(Tmx_clim,[251,1]);
%%T_glb_ssp126 = [tas_historical;tas_ssp126];
%%T_glb_ssp126 = T_glb_ssp126 - repmat(Tmx_clim,[251,1]);
%t_est = 1:(365*251);
%coeff_rom = rom_LE_glbT_lin_mean_var_pred(coeff_ssp(:,:,mem_start:mem_end),T_glb,T_glb_ssp119,t_est,t_est);
%clear coeff_ssp;
disp('start loading coeff_rom');
coeff_rom_full = zeros(l1,n_eig);
for k = 1:(n_eig)/50
    disp(k);
    istart = (k-1)*50+1;
    iend = k*50;
    if(scenario == 'ssp585')
        filename = strcat(folder_coeff,file_coeff,scenario,'_ens_',num2str(n_ens,'%02d'),'_',num2str(k,'%02d'),'.mat');
        load(filename,'coeff_rom');
        coeff_rom_full(:,istart:iend) = coeff_rom(:,:,n_seed);
    else
        filename = strcat(folder_coeff,file_coeff,'ssp585','_ens_',num2str(n_ens,'%02d'),'_',num2str(k,'%02d'),'.mat');
        load(filename,'coeff_rom');
        daystart = (1850 - 1850)*365 + 1;
        dayend   = (2014 - 1850 + 1)*365;
        coeff_rom_full(daystart:dayend,istart:iend) = coeff_rom(daystart:dayend,:,n_seed);

        filename = strcat(folder_coeff,file_coeff,scenario,'_ens_',num2str(n_ens,'%02d'),'_',num2str(k,'%02d'),'.mat');
        load(filename,'coeff_rom');
        daystart = (2015 - 1850)*365 + 1;
        dayend   = (2100 - 1850 + 1)*365;
        coeff_rom_full(daystart:dayend,istart:iend) = coeff_rom(:,:,n_seed);
    end
end
coeff_rom = coeff_rom_full;
disp('finish loading coeff_rom');

%% load PCA modes
inname = strcat('pca_historical/modes_Tmx_365skip5_global_historical_ens.mat');
load(inname,'w_hat','lambda');
mode = w_hat;
clear w_hat;

%% post processing
for statsYear = startYear:skipYear:endYear
for l = n_seed
member = strcat('r',num2str(l,'%d'),'i1p1f1');
disp(l);
    %k = statsYear - startYear + 1;
    currentDate = datetime(statsYear,1,1);

    % collect samples from the local neighborhood
    Tmx_c = zeros(n_lonc,n_latc,92*lat_skip*lon_skip*skipYear);
    nn = 0;
    nlength = 92*lat_skip*lon_skip;

    for readYear = statsYear:statsYear+skipYear-1    
        disp(readYear);
        % PCA reconstruction
        Tmx_est = zeros(n_lon,n_lat,365);
        day_start = 152;
        day_end = 243;
        %if(readYear < 2015)
        k = readYear - 1850 + 1;
        tic;
        for j=day_start:day_end
        for i_mode = 1:n_eig
            indx = 365*(k-1)+j;
            Tmx_est(:,:,j) = Tmx_est(:,:,j) + coeff_rom(indx,i_mode) * mode(:,:,i_mode);
        end
        end
        % add climatological mean
        Tmx_est(:,:,day_start:day_end) = Tmx_est(:,:,day_start:day_end) + ...
                                     Tmx_mean(:,:,day_start:day_end);
        Tmx = Tmx_est(:,:,day_start:day_end);
        clear Tmx_est;
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
        
        nn = nn + 1;
     end

     % mean and standard deviation
     TMX_avg = mean(Tmx_c,3);
     TMX_std = std(Tmx_c,0,3);
     
     if(statsYear + skipYear - 1 < 2015)
         outname = strcat(folder_output,'stats_',file_output,'_10year_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     else
         outname = strcat(folder_output,'stats_',file_output,'_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     end
     save(outname,'TMX_avg','TMX_std','lonc','latc','skipYear');

     %% compute quantiles
     n_samp = 92*lat_skip*lon_skip*skipYear;
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
         outname = strcat(folder_output,'qntl_',file_output,'_10year_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     else
         outname = strcat(folder_output,'qntl_',file_output,'_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     end
     save(outname,'qntl_75','qntl_90','qntl_95','qntl_975');

     %% extract local PDF
     n_samp = 92*lat_skip*lon_skip*skipYear;
     Tmx_local = zeros(n_samp,19);
     city_coords = readtable('stats_10year/city_coords.csv');
     for i=1:19
         if(city_coords.lon(i) < 0)
             lon_tmp = city_coords.lon(i) + 360;
         else
             lon_tmp = city_coords.lon(i);
         end
         if(lon_tmp > 359.5)
             lon_tmp = 0;
         end
         idx_lon = find(abs(lonc - lon_tmp) < 0.9375*lon_skip);
         idx_lat = find(abs(latc - city_coords.lat(i)) < 0.9325*lat_skip);
         disp([lonc(idx_lon),latc(idx_lat)]);
         tmp = reshape(Tmx_c(idx_lon,idx_lat,:),[n_samp,1]);
         Tmx_local(:,i) = tmp; 
     end

     if(statsYear + skipYear - 1 < 2015)
         outname = strcat(folder_output,'localPDF_',file_output,'_10year_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     else
         outname = strcat(folder_output,'localPDF_',file_output,'_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     end
     save(outname,'Tmx_local','city_coords');
end
end
