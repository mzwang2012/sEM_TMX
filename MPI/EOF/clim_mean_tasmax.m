clear;
close all;

%% Initialization
flag = 0;   % 0 for computing climatological mean
            % 1 for computing global averaged variance

%folder = '/net/fs09/d0/lutjens/raw/CMIP6/MPI-ESM1-2-LR/r1i1p1f1/historical/tasmax/250_km/day/';
folder = '/net/fs09/d0/lutjens/raw/CMIP6/MPI-ESM1-2-LR/';

n_ens = 50;

scenario = '';
variable = 'tasmax';

startYear = 1850;
endYear = 2014;
skipDay = 5;
n_years = endYear - startYear + 1;

% load coordinates
filename = strcat(folder,'r1i1p1f1/historical/tasmax/250_km/day/1850/CMIP6_MPI-ESM1-2-LR_r1i1p1f1_historical_tasmax_250_km_day_gn_1850.nc');
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
n_lat = length(lat);
n_lon = length(lon);

if(flag == 0)
    %% compute climatological mean
    Tmx_mean = zeros(n_lon,n_lat,365);
    for l = 1:n_ens
        member = strcat('r',num2str(l,'%d'),'i1p1f1');
        disp(member);
        for k = 1:n_years
            currentYear = startYear + k - 1;
            disp(currentYear);
            % read data from CMIP6, returns Tmx(n_lon,n_lat,365)
            Tmx = MPI_LR_read_tasmax(folder,member,scenario,variable,currentYear,n_lon,n_lat);
        
            Tmx_mean = Tmx_mean + Tmx;
        end
    end
    
    Tmx_mean = Tmx_mean / n_years / n_ens;
    
    save('pca_global/Tmx_clim_historical.mat','Tmx_mean','startYear','endYear','scenario','variable','lon','lat');

end
