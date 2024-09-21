%% load tasmax data from cmip6 MPI-LR .nc file
function Tmx = MPI_LR_read_tasmax(folder,member,scenario,variable,currentYear,n_lon,n_lat)

    % differentiate historical and future scenarios
    if(currentYear < 1850 || currentYear > 2100)
        disp('Data before 1850 or after 2100 are not available!!!');
        return;
    elseif(currentYear >= 1850 && currentYear < 2015)
        filename = strcat(folder,member,'/historical/',variable,...
            '/250_km/day/',num2str(currentYear),'/CMIP6_MPI-ESM1-2-LR_',...
            member,'_historical_',variable,'_250_km_day_gn_',num2str(currentYear),'.nc');
    elseif(currentYear >= 2015 && currentYear <= 2100)
        filename = strcat(folder,member,'/',scenario,'/',variable,...
            '/250_km/day/',num2str(currentYear),'/CMIP6_MPI-ESM1-2-LR_',...
            member,'_',scenario,'_',variable,'_250_km_day_gn_',num2str(currentYear),'.nc');
    end

    % load Tmx data
    if(leapyear(currentYear))
      disp(['leap year at ',num2str(currentYear)]);

      start = [1 1 1];
      count = [n_lon,n_lat,366];
      
      Tmx = ncread(filename,variable,start,count);
      %disp(['shape of Tmx: ',num2str(size(Tmx))]);

      % remove 02/29 data
      Tmx(:,:,60) = [];
      %disp(['shape of Tmx: ',num2str(size(Tmx))]);
   
   else
      start = [1 1 1];
      count = [n_lon,n_lat,365];
      
      Tmx = ncread(filename,variable,start,count);
      %disp(['shape of Tmx: ',num2str(size(Tmx))]);
   end
