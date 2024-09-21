clear;
close all;

variable = 'tas';

%% load parameters from nc file
%folder = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/historical/';
%filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_18500101-18991231.nc');
%filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_19000101-19491231.nc');
%filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_19500101-19991231.nc');
%filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_20000101-20141231.nc');
%folder = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/ssp126/';
%filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_ssp126_r1i1p1f2_gr_20150101-20641231.nc');
%filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_ssp126_r1i1p1f2_gr_20650101-21001231.nc');
folder = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/ssp585/';
%filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_ssp585_r1i1p1f2_gr_20150101-20641231.nc');
filename = strcat(folder,variable,'/tas_day_CNRM-CM6-1-HR_ssp585_r1i1p1f2_gr_20650101-21001231.nc');
lat = ncread(filename,'lat');
n_lat = length(lat);
lon = ncread(filename,'lon');
n_lon = length(lon);
height = ncread(filename,'height');
time = ncread(filename,'time');
% time_bounds = ncread(filename,'time_bounds');

%% degree to radian
theta = 90 - lat;
theta_rad = theta / 180 * pi;
phi_rad = lon / 180 *pi;
Mat_sin = reshape(sin(theta_rad),[1,n_lat]);

%% Days to date
% Assuming time is your vector storing the # of days since 01/01/1850
refDate = datetime(1850,1,1);
dates = refDate + days(time-0.5); % Convert to datetime array
startDate = refDate + days(min(time)); 
endDate = refDate + days(max(time)); % Find the end date based on the maximum value in time

% Determine the range of years
startYear = year(startDate);
endYear = year(endDate);

% Extract year, month, and day
%years = year(dates);
%months = month(dates);
%days = day(dates);

n_years = endYear - startYear + 1;
disp(n_years);

% skip 02/29 in leap years
T_glb_mean = zeros(n_years*365,1);

%% loop over years
istart = 0;
iend = 0;
for k = 1:n_years
   currentYear = startYear + k - 1;
   disp(currentYear);
   istart = iend + 1;
   iend = iend + 365;

   if(leapyear(currentYear))
      disp(['leap year at ',num2str(currentYear)]);

      start = [1 1 istart];
      count = [n_lon,n_lat,366];
      
      Tmx = ncread(filename,variable,start,count);
      disp(['shape of Tmx: ',num2str(size(Tmx))]);

      % remove 02/29 data
      Tmx(:,:,60) = [];
      disp(['shape of Tmx: ',num2str(size(Tmx))]);
   
   else
      start = [1 1 istart];
      count = [n_lon,n_lat,365];
      
      Tmx = ncread(filename,variable,start,count);
      disp(['shape of Tmx: ',num2str(size(Tmx))]);
   end

   %% global mean: 2D integration \int 
   Tmx_lat = trapz(phi_rad,Tmx,1);
   Tmx_mean = trapz(theta_rad,Tmx_lat.*Mat_sin,2) / (-4*pi);
   
   T_glb_mean(istart:iend) = squeeze(Tmx_mean);
end

%save('global_stats/historical_tas_mean_4.mat','time','startYear','endYear','T_glb_mean');
save('global_stats/ssp585_tas_mean_2.mat','time','startYear','endYear','T_glb_mean');
