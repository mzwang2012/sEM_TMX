%% load Tmx data from cmip6 .nc files
function Tmx = cmip6_read_Tmx(folder,scenario,variable,currentYear,n_lon,n_lat)
    currentDate = datetime(currentYear,1,1);
    % differentiate which file to load
    if(currentYear < 1850)
        disp('There is a bug in the code!!!');
        return;
    elseif(currentYear >= 1850 && currentYear < 1900)
        filename = strcat(folder,'historical/',variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_18500101-18991231.nc');
        refDate = datetime(1850,1,1);        
    elseif(currentYear >= 1900 && currentYear < 1950)
        filename = strcat(folder,'historical/',variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_19000101-19491231.nc');
        refDate = datetime(1900,1,1); 
    elseif(currentYear >= 1950 && currentYear < 2000)
        filename = strcat(folder,'historical/',variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_19500101-19991231.nc');
        refDate = datetime(1950,1,1); 
    elseif(currentYear >= 2000 && currentYear < 2015)
        filename = strcat(folder,'historical/',variable,...
            '/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_20000101-20141231.nc');
        refDate = datetime(2000,1,1); 
    elseif(currentYear >= 2015 && currentYear < 2065)
        filename = strcat(folder,scenario,'/',variable,...
            '/tasmax_day_CNRM-CM6-1-HR_',scenario,'_r1i1p1f2_gr_20150101-20641231.nc');
        refDate = datetime(2015,1,1);        
    elseif(currentYear >= 2065 && currentYear < 2101)
        filename = strcat(folder,scenario,'/',variable,...
            '/tasmax_day_CNRM-CM6-1-HR_',scenario,'_r1i1p1f2_gr_20650101-21001231.nc');
        refDate = datetime(2065,1,1); 
    end

    % count the range of indices to load
    duration = days(currentDate - refDate);
    istart = duration + 1;
    %disp(istart);

    % load Tmx data
    if(leapyear(currentYear))
      disp(['leap year at ',num2str(currentYear)]);

      start = [1 1 istart];
      count = [n_lon,n_lat,366];
      
      Tmx = ncread(filename,variable,start,count);
      %disp(['shape of Tmx: ',num2str(size(Tmx))]);

      % remove 02/29 data
      Tmx(:,:,60) = [];
      %disp(['shape of Tmx: ',num2str(size(Tmx))]);
   
   else
      start = [1 1 istart];
      count = [n_lon,n_lat,365];
      
      Tmx = ncread(filename,variable,start,count);
      %disp(['shape of Tmx: ',num2str(size(Tmx))]);
   end
