clear;
close all;

% get the name of the currently running script
currentFile = mfilename;
fileIdx = str2num(currentFile(end-1:end));
disp(fileIdx);

%% path to CMIP6 data
folder = '/net/fs06/d3/lutjens/bc3/data/raw/CMIP6/MPI-ESM1-2-LR/';
scenario = 'ssp119';
variable = 'tasmax';
n_variables = 1;    % # of variables

mem_start = 1; 
mem_end = 50;
n_ens = (mem_end - mem_start)+1;

startYear = 1850;
endYear = 2100;
n_years = endYear - startYear + 1;
n_days =  n_years * 365;

% load coordinates
load('pca_historical/Tmx_clim_historical.mat','lon','lat','Tmx_mean');
n_lat = length(lat);
n_lon = length(lon);

%% load climatological mean and PCA modes
inname = strcat('pca_historical/modes_Tmx_365skip5_global_historical_ens.mat');
load(inname,'w_hat','lambda');
mode = w_hat;
clear w_hat;
offset = 0;

[~,~,n_eig] = size(mode);

% range of mode indices
startMode = (fileIdx-offset-1)*50+1;
endMode = (fileIdx-offset)*50;

coeff = zeros(n_days,endMode - startMode + 1,n_ens);

S_vec = sqrt(sin((90 - lat)/180*pi));
S_mat = repmat(reshape(S_vec,[1,n_lat]),[n_lon,1]);
%% projection CMIP6 data onto PCA modes
for k = 1:n_years
for l = 1:n_ens
    member = strcat('r',num2str(l,'%d'),'i1p1f1');
    disp(member);
    currentYear = startYear + k - 1;
    disp(currentYear);
    currentDate = datetime(currentYear,1,1);
     
    Tmx = MPI_LR_read_tasmax(folder,member,scenario,variable,currentYear,n_lon,n_lat);

   % remove climatological mean
    Tmx = Tmx - Tmx_mean;

   % remove cproject onto PCA modes
   for j = startMode:endMode
       jndx = j - startMode + 1;
       for i = 1:365
           indx = 365*(k-1)+i;
           proj_Tmx =  sum(sum(Tmx(:,:,i) .* mode(:,:,j) .* S_mat(:,:,1).^2,1),2);
           coeff(indx,jndx,l) = proj_Tmx;
       end
   end

end
end
outname = strcat('pca_historical/Tmx_365_global_',scenario,'_coeff_',num2str(fileIdx,'%02d'),'.mat');
save(outname,'startYear','endYear','startMode','endMode','coeff','mem_start','mem_end');
