clear;
close all;

%currentFile = mfilename;
%fileIdx = str2num(currentFile(end-2:end));
%disp(fileIdx);
%
%mem_id = fileIdx;

%% path to CMIP6 data
folder = '/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/';
scenario = 'ssp126';
variable = 'tasmax';

mem_start = 1; 
mem_end = 50;
n_ens = (mem_end - mem_start)+1;

lat_skip = 1;
lon_skip = 1;
%% load climatological mean and PCA modes
load('pca_historical/Tmx_clim_historical.mat','lon','lat');
n_lat = length(lat);
n_lon = length(lon);

% load('pca_ERA5/Tmx_clim_1979_1998.mat','Tmx_mean');

startYear = 2010; %1850;
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

%% post processing
for statsYear = startYear:skipYear:endYear
for l = mem_start:mem_end
member = strcat('r',num2str(l,'%d'),'i1p1f1');
disp(member);
    %k = statsYear - startYear + 1;
    currentDate = datetime(statsYear,1,1);

    % collect samples from the local neighborhood
    Tmx_c = zeros(n_lonc,n_latc,92*lat_skip*lon_skip*skipYear);
    nn = 0;
    nlength = 92*lat_skip*lon_skip;

    for readYear = statsYear:statsYear+skipYear-1    
        disp(readYear);
        % load true field
        Tmx = MPI_LR_read_tasmax(folder,member,scenario,variable,readYear,n_lon,n_lat);
        Tmx = Tmx(:,:,152:243);
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
         outname = strcat('stats_10year/stats_10year_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     else
         outname = strcat('stats_10year/stats_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
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
         outname = strcat('stats_10year/qntl_10year_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     else
         outname = strcat('stats_10year/qntl_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
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
         outname = strcat('stats_10year/localPDF_10year_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     else
         outname = strcat('stats_10year/localPDF_10year_',scenario,'_',num2str(statsYear,'%d'),'_',num2str(l,'%03d'),'.mat');
     end
     save(outname,'Tmx_local','city_coords');
end
end
