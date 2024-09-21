%% project TMX data onto EOFs
clear;
close all;

%% path to CMIP6 data
% folder = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/historical/';
folder = 'download/';
variable = 'tasmax';

%% load climatological mean and PCA modes
load('pca_ERA5/Tmx_365_global_1979_1998_modes.mat');
[n_lon,n_lat,n_eig] = size(mode);

load('pca_ERA5/Tmx_clim_1979_1998.mat','Tmx_mean');

startYear = 1850;
endYear = 1850; %2014;
n_years = endYear - startYear + 1;
n_days =  n_years * 365;

% range of mode indices
startMode = 1;
endMode = 1;

coeff = zeros(n_days,endMode - startMode + 1);

S_vec = sqrt(sin((90 - lat)/180*pi));
S_mat = repmat(reshape(S_vec,[1,n_lat]),[n_lon,1]);
%% projection CMIP6 data onto PCA modes
for k = 1:n_years
    currentYear = startYear + k - 1;
    disp(currentYear);
    currentDate = datetime(currentYear,1,1);

    % differentiate which file to load
    if(currentYear < 1850)
        disp('There is a bug in the code!!!');
        return;
    elseif(currentYear >= 1850 && currentYear < 1900)
        filename = strcat(folder,variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_18500101-18991231.nc');
        refDate = datetime(1850,1,1);        
    elseif(currentYear >= 1900 && currentYear < 1950)
        filename = strcat(folder,variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_19000101-19491231.nc');
        refDate = datetime(1900,1,1); 
    elseif(currentYear >= 1950 && currentYear < 2000)
        filename = strcat(folder,variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_19500101-19991231.nc');
        refDate = datetime(1950,1,1); 
    elseif(currentYear >= 2000 && currentYear < 2015)
        filename = strcat(folder,variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_20000101-20141231.nc');
        refDate = datetime(2000,1,1); 
    end

    % count the range of indices to load
    duration = days(currentDate - refDate);
    istart = duration + 1;
    disp(istart);

    % load Tmx data
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

   % remove climatological mean and project onto PCA modes
   Tmx = Tmx - Tmx_mean;
   for j = startMode:endMode
       jndx = j - startMode + 1;
       for i = 1:365
           indx = 365*(k-1)+i;
           coeff(indx,jndx) = sum(sum(Tmx(:,:,i) .* mode(:,:,j) .* S_mat(:,:,1).^2,1),2);
       end
   end

end
