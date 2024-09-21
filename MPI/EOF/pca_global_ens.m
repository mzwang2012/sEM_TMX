clear;
close all;

%currentFile = mfilename;
%fileIdx = str2num(currentFile(end-2:end));
%disp(fileIdx);
%
%mem_id = fileIdx;

%% Initialization
flag = 3;   % 1 for computing covariance matrix
            % 2 for averaging xxt matrix over large ensembles
            % 3 for eigen decomposition

% path to CMIP6 data
% folder = 'mnt/';
folder = '/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/';
scenario = 'historical';
variable = 'tasmax';

%mem_start = mem_id; 
%mem_end = mem_id;
mem_start = 1; 
mem_end = 50;
n_ens = (mem_end - mem_start)+1;

startYear = 1850;
endYear = 2014;
skipDay = 5;
n_eig = 2000;
n_years = endYear - startYear + 1;
%n_days =  n_years * 365;

% load coordinates
load('pca_historical/Tmx_clim_historical.mat','lon','lat');
n_lat = length(lat);
n_lon = length(lon);

if(flag == 1)
    %% compute covariance matrix for PCA
    % initialization
    n_days = length(1:skipDay:365);  % # of days to be sampled every year
    n_snap = n_years*n_days;
    S_vec = sqrt(sin((90 - lat)/180*pi));
    S_mat = repmat(reshape(S_vec,[1,n_lat]),[n_lon,1]);
    nlength = n_lon * n_lat;    
    mat_xxt = zeros(nlength,nlength);

    load('pca_historical/Tmx_clim_historical.mat','Tmx_mean');

    istart = 0;
    iend = 0;
    % restart
    % inname = strcat('pca_historical/matxtx_Tmx_365skip2_global_historical.mat');
    % load(inname);
    % tmpYear = year1;  
    for l = mem_start:mem_end
        member = strcat('r',num2str(l,'%d'),'i1p1f1');
        disp(member);
        for year1 = startYear:endYear
            disp(year1);
            Tmx = MPI_LR_read_tasmax(folder,member,scenario,variable,year1,n_lon,n_lat);
            Tmx = (Tmx - Tmx_mean) .* S_mat;
            Tmx_1 = reshape(Tmx(:,:,1:skipDay:365),[nlength,n_days]);
            clear Tmx;
            
            mat_xxt = mat_xxt + Tmx_1 * Tmx_1';
        end
    end
  
    mat_xxt = mat_xxt / n_ens / n_snap;
    outname = strcat('pca_historical/matxxt_Tmx_365skip5_global_historical_mem_',num2str(mem_id,'%03d'),'.mat');
    save(outname,'mat_xxt','startYear','endYear','variable','scenario','skipDay',...
         'mem_start','mem_end','-v7.3');
elseif(flag == 2)
    nlength = n_lon * n_lat;    
    mat_xxt_ens = zeros(nlength,nlength);
    for mem_id = mem_start:mem_end
        disp(mem_id);
        inname = strcat('pca_historical/matxxt_Tmx_365skip5_global_historical_mem_',num2str(mem_id,'%03d'),'.mat');
        load(inname,'mat_xxt');
        mat_xxt_ens = mat_xxt_ens + mat_xxt;
    end
    outname = strcat('pca_historical/matxxt_Tmx_365skip5_global_historical_ens.mat');
    save(outname,'mat_xxt','startYear','endYear','variable','scenario','skipDay',...
         'mem_start','mem_end','-v7.3');

elseif(flag == 3)

    %% compute PCA modes
    % initialization
    n_days = length(1:skipDay:365);  % # of days to be sampled every year
    n_snap = n_years*n_days;
    S_vec = sqrt(sin((90 - lat)/180*pi));
    S_mat = repmat(reshape(S_vec,[1,n_lat]),[n_lon,1]);
    nlength = n_lon * n_lat;    
    inname = strcat('pca_historical/matxxt_Tmx_365skip5_global_historical_ens.mat');
    load(inname,'mat_xxt');
    disp('start eig');
    [V,D] = eig(mat_xxt);
    disp('end eig');
    [lambda,ind] = sort(diag(abs(D)),'descend');
    % Note: there should be 365 zero eigenvalues because we are subtracting climatological mean
    V = V(:,ind);

    w_hat = reshape(V(:,1:n_eig),[n_lon,n_lat,n_eig]);
    for i=1:n_eig
        w_hat(:,:,i) = w_hat(:,:,i) ./ S_mat(:,:,1);
    end

    % normalization
    disp('start normalization');
    for i=1:n_eig
        tmp = w_hat(:,:,i);
        const = sum(sum((tmp .* S_mat(:,:,1)).^2,1),2);
        const = sqrt(const);
        disp(const);
        w_hat(:,:,i) = w_hat(:,:,i) / const;
    end
    %if(n_eig > 500)
    %    w_hat_2 = w_hat(:,:,501:n_eig);
    %    w_hat = w_hat(:,:,1:500);
    %end
    %outname = strcat('pca_historical/modes_Tmx_365skip5_global_historical.mat');
    %save(outname,'lambda','w_hat','w_hat_2','n_eig','variable','scenario','startYear','endYear','skipDay','-v7.3');
    outname = strcat('pca_historical/modes_Tmx_365skip5_global_historical_ens.mat');
    save(outname,'lambda','w_hat','n_eig','variable','scenario','startYear','endYear','skipDay',...
                 'mem_start','mem_end','-v7.3');
end

