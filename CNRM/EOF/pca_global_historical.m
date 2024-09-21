%% compute EOFs / PCA modes

clear;
close all;

%% Initialization
flag = 2;   % 0 for computing climatological mean
            % 1 for computing covariance matrix
            % 2 for eigen decomposition

% path to CMIP6 data
folder = '/net/fs06/d3/CMIP6/CNRM-CM6-1-HR/';
scenario = 'historical';
variable = 'tasmax';

startYear = 1850;
endYear = 2014;
skipDay = 5;
n_eig = 2000;
n_years = endYear - startYear + 1;
%n_days =  n_years * 365;

% load coordinates
filename = strcat(folder,'historical/tasmax/tasmax_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_18500101-18991231.nc');
lat = ncread(filename,'lat');
n_lat = length(lat);
lon = ncread(filename,'lon');
n_lon = length(lon);
height = ncread(filename,'height');

if(flag == 0)
    %% compute climatological mean
    Tmx_mean = zeros(n_lon,n_lat,365);
    for k = 1:n_years
        currentYear = startYear + k - 1;
        disp(currentYear);
        % read data from CMIP6, returns Tmx(n_lon,n_lat,365)
        Tmx = cmip6_read_Tmx(folder,scenario,variable,currentYear,n_lon,n_lat);
    
        Tmx_mean = Tmx_mean + Tmx;
    end
    
    Tmx_mean = Tmx_mean / n_years;
    
    save('pca_historical/Tmx_clim_historical.mat','Tmx_mean','startYear','endYear','scenario','variable');
    
elseif(flag == 1)
    %% compute covariance matrix for PCA
    % initialization
    n_days = length(1:skipDay:365);  % # of days to be sampled every year
    n_snap = n_years*n_days;
    S_vec = sqrt(sin((90 - lat)/180*pi));
    S_mat = repmat(reshape(S_vec,[1,n_lat]),[n_lon,1]);
    nlength = n_lon * n_lat;    
    mat_xtx = zeros(n_snap,n_snap);

    load('pca_historical/Tmx_clim_historical.mat','Tmx_mean');

    istart = 0;
    iend = 0;
    % restart
    % inname = strcat('pca_historical/matxtx_Tmx_365skip2_global_historical.mat');
    % load(inname);
    % tmpYear = year1;    
    for year1 = startYear:endYear
    % for year1 = tmpYear:endYear
        disp(year1);
        Tmx = cmip6_read_Tmx(folder,scenario,variable,year1,n_lon,n_lat);
        
        % subtract the mean field
        % 2D matrix S_mat will be multiplied by each slice of 3D Tmx
        Tmx = (Tmx - Tmx_mean) .* S_mat;
        
        Tmx_1 = reshape(Tmx(:,:,1:skipDay:365),[nlength,n_days]);
        clear Tmx;
        
        jstart = istart;
        jend = iend;
        
        istart = iend + 1;
        iend = iend + n_days;
        
        for year2 = year1:endYear
            disp([year1,year2]);
            Tmx = cmip6_read_Tmx(folder,scenario,variable,year2,n_lon,n_lat);
            
            % subtract the mean field
            Tmx = (Tmx - Tmx_mean) .* S_mat;
            Tmx_2 = reshape(Tmx(:,:,1:skipDay:365),[nlength,n_days]);
            clear Tmx;
%             istart = (year  -year_start  )*ndays_1+1;
%             iend   = (year  -year_start+1)*ndays_1;
%             jstart = (year2 -year_start  )*ndays_2+1;
%             jend   = (year2 -year_start+1)*ndays_2;
            
            jstart = jend + 1;
            jend = jend + n_days;
            
            mat_xtx(istart:iend,jstart:jend) = Tmx_1' * Tmx_2;
            mat_xtx(jstart:jend,istart:iend) = Tmx_2' * Tmx_1;
            
        end
        outname = strcat('pca_historical/matxtx_Tmx_365skip2_global_historical.mat');
        save(outname,'mat_xtx','startYear','endYear','variable','scenario','skipDay',...
             'istart','iend','jstart','jend','year1','year2','-v7.3');
    end
    outname = strcat('pca_historical/matxtx_Tmx_365skip2_global_historical.mat');
    save(outname,'mat_xtx','startYear','endYear','variable','scenario','skipDay','-v7.3');

elseif(flag == 2)
    %% compute PCA modes
    % initialization
    n_days = length(1:skipDay:365);  % # of days to be sampled every year
    n_snap = n_years*n_days;
    S_vec = sqrt(sin((90 - lat)/180*pi));
    S_mat = repmat(reshape(S_vec,[1,n_lat]),[n_lon,1]);
    nlength = n_lon * n_lat;    
    inname = strcat('pca_historical/matxtx_Tmx_365skip5_global_historical.mat');
    load(inname,'mat_xtx');
    disp('start eig');
    [V,D] = eig(mat_xtx);
    disp('end eig');
    [lambda,ind] = sort(diag(abs(D)),'descend');
    % Note: there should be 365 zero eigenvalues because we are subtracting climatological mean
    V = V(:,ind);

    load('pca_historical/Tmx_clim_historical.mat','Tmx_mean');

    istart = 0;
    iend = 0;
    w_hat = zeros(nlength,n_eig);
    for year1 = startYear:endYear
        disp(year1);
        Tmx = cmip6_read_Tmx(folder,scenario,variable,year1,n_lon,n_lat);
        
        % subtract the mean field
        % 2D matrix S_mat will be multiplied by each slice of 3D Tmx
        Tmx = (Tmx - Tmx_mean) .* S_mat;
        
        Tmx_1 = reshape(Tmx(:,:,1:skipDay:365),[nlength,n_days]);
        clear Tmx;
        
        
        istart = iend + 1;
        iend = iend + n_days;
        w_hat = w_hat + Tmx_1 * V(istart:iend,1:n_eig);
    end

    w_hat = w_hat * diag((1./sqrt(lambda(1:n_eig))));
    w_hat = reshape(w_hat,[n_lon,n_lat,n_eig]);
    for i=1:n_eig
        w_hat(:,:,i) = w_hat(:,:,i) ./ S_mat(:,:,1);
    end

    % normalization
    disp('start normalization');
    for i=1:n_eig
        tmp = w_hat(:,:,i);
        const = sum(sum((tmp .* S_mat(:,:,1)).^2,1),2);
        const = sqrt(const);
        w_hat(:,:,i) = w_hat(:,:,i) / const;
    end
    
    w_hat_3 = w_hat(:,:,1001:1500);
    w_hat_4 = w_hat(:,:,1501:n_eig);
    outname = strcat('pca_historical/modes2_Tmx_365skip5_global_historical.mat');
    save(outname,'lambda','w_hat_3','w_hat_4','n_eig','variable','scenario','startYear','endYear','skipDay','-v7.3');

    w_hat_2 = w_hat(:,:,501:1000);
    w_hat = w_hat(:,:,1:500);
    outname = strcat('pca_historical/modes_Tmx_365skip5_global_historical.mat');
    save(outname,'lambda','w_hat','w_hat_2','n_eig','variable','scenario','startYear','endYear','skipDay','-v7.3');
end

