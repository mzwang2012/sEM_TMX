clear;
close all;

%% Initialization
flag = 0;   % 0 for computing global mean
            % 1 for computing global averaged variance

%folder = '/net/fs09/d0/lutjens/raw/CMIP6/MPI-ESM1-2-LR/r1i1p1f1/historical/tasmax/250_km/day/';
folder = '/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/';

n_ens = 50;

scenario = '';
variable = 'tas';

startYear = 1850;
endYear = 2014;
skipDay = 1;
n_years = endYear - startYear + 1;

% load coordinates
load('pca_historical/Tmx_clim_historical.mat','lon','lat');
n_lat = length(lat);
n_lon = length(lon);
%% degree to radian
theta = 90 - lat;
theta_rad = theta / 180 * pi;
phi_rad = lon / 180 *pi;
Mat_sin = reshape(sin(theta_rad),[1,n_lat]);

if(flag == 0)
    %% compute global mean
    T_glb_mean = zeros(n_years*365,1);
    for l = 1:n_ens
        member = strcat('r',num2str(l,'%d'),'i1p1f1');
        disp(member);
        istart = 0;
        iend = 0;
        T_glb_mean_tmp = zeros(n_years*365,1);
        for k = 1:n_years
            istart = iend + 1;
            iend = iend + 365;
            currentYear = startYear + k - 1;
            disp(currentYear);
            % read data from CMIP6, returns Tmx(n_lon,n_lat,365)
            Tmx = MPI_LR_read_tasmax(folder,member,scenario,variable,currentYear,n_lon,n_lat);
        
            %% global mean: 2D integration \int 
            Tmx_lat = trapz(phi_rad,Tmx,1);
            Tmx_mean = trapz(theta_rad,Tmx_lat.*Mat_sin,2) / (-4*pi);
            disp(Tmx_mean(1));
            T_glb_mean_tmp(istart:iend) = squeeze(Tmx_mean);
        end
        T_glb_mean = T_glb_mean + T_glb_mean_tmp;
    end
    
    Tmx_glb_mean = Tmx_glb_mean / n_ens;
    
    save('global_stats/historical_tas_mean.mat','time','startYear','endYear','T_glb_mean');

end
