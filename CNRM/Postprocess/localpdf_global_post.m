clear;
close all;

choose_data = 0; % 1 for loading data from mounted folder; 
choose_post = 1; % 1 for local PDF, averaged in 4 years;
                 % 2 for local PDF, averaged in 20 years;
                 % 3 for local PDF, at mode maxima


year_skip = 10;  % number of years sampled altogether
flag_print = 0;

lat_skip = 2;
lon_skip = 2;
n_eig = 50;

%% load grid and get the coarse grid
folder = 'download/';
variable = 'tasmax';

filename = 'pca_ERA5/Tmx_365_global_1979_1998_modes.mat';
load(filename,'lat','lon');

n_lat = length(lat);
n_lon = length(lon);

n_latc = floor(n_lat/lat_skip);
n_lonc = floor(n_lon/lon_skip);

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
if(choose_post == 1)
load coastlines;
% min_x_store = zeros(19,1);
% max_x_store = zeros(19,1);
% max_y_store = zeros(19,1);
% load('pdf_axis_range.mat');

for i= 11%:19 %1:4

% min_x = 10;
% max_x = 35;
% max_y = 0.12;

% scenario = '';
% for year = 2000 %1850:year_skip:2000 %2010:year_skip:2090
scenario = 'ssp585_';
for year = 2090 %2010:year_skip:2090 
    filename = strcat('stats_autumn/localPDF_10year_',scenario,num2str(year,'%d'),'.mat');
    load(filename);
    tmp = Tmx_local(:,i) - 273.15;
    data_true = tmp;
%     filename = strcat('stats_10year/localPDF_est_10year_',num2str(year,'%d'),'.mat');
%     filename = strcat('stats_10year/rom_glbT_true_mean_var/',...
%         'localPDF_rom_10year_',scenario,num2str(year,'%d'),'_001.mat');
%     load(filename);
%     tmp = Tmx_local(:,i) - 273.15;
%     data_est = tmp;
    in_buffer = zeros(length(tmp),10);
    for kk=1:10
        filename = strcat('stats_autumn/localPDF_rom_10year_',...
            scenario,num2str(year,'%d'),...,
            '_',num2str(kk,'%03d'),'.mat');

%         filename = strcat('stats_noice/localPDF_rom_noice_',...
%             scenario,num2str(year,'%d'),...,
%             '_',num2str(kk,'%03d'),'.mat');
        load(filename);
        in_buffer(:,kk) = Tmx_local(:,i);
    end
    data_rom = in_buffer - 273.15;
    
    [f,xi] = ksdensity(data_true);
    epsilon = 1e-8;
    nboot = 10;
    pdf_boot = zeros(nboot,numel(xi));
    for k=1:nboot
        sample = datasample(data_true,numel(data_true));
        pdf_boot(k,:) = ksdensity(sample,xi);
    end
    pdf_error = std(pdf_boot);
    f_low = f - pdf_error;
    f_upp = f + pdf_error;
    f_low(f_low < 0) = epsilon;
    
%     [f_est,xi_est] = ksdensity(data_est);
    [f_rom,xi_rom] = ksdensity(data_rom(:));
% 
%     % PDF of ROM 
    f_rom_mrlz = zeros(numel(xi_rom),10);
    for l=1:10
        f_rom_mrlz(:,l) = ksdensity(data_rom(:,l),xi_rom);
    end
    f_rom_error = std(f_rom_mrlz,[],2)';
    f_rom_low = f_rom - f_rom_error;
    f_rom_upp = f_rom + f_rom_error;
    f_rom_low(f_rom_low < 0) = epsilon;
    mu_true = mean(data_true);
%     mu_est = mean(data_est);
    mu_rom = mean(data_rom(:));

    figure;
    fill([xi,fliplr(xi)],[f_low,fliplr(f_upp)],...
        [1 0.7 0.7],'EdgeColor','none');
    hold on;
    fill([xi_rom,fliplr(xi_rom)],[f_rom_low,fliplr(f_rom_upp)],...
        [0.7 0.7 1],'EdgeColor','none');

    plot(xi,f,'r--','linewidth',2);
    hold on;  
%     plot(xi_est,f_est,'k','linewidth',2);
    plot(xi_rom,f_rom,'b','linewidth',2);

    plot([mu_true,mu_true],[0,max(f)],'r:','linewidth',2);
%     plot([mu_est,mu_est],[0,max(f_est)],'k:','linewidth',2);
    plot([mu_rom,mu_rom],[0,max(f_rom)],'b:','linewidth',2);
    %     plot(xi_tm,f_tm,'b','linewidth',2);
%     if(year == 1850)
        min_x = min([min(xi_rom),min(xi)]);
        max_x = max([max(xi_rom),max(xi)])+2;
        max_y = max([max(f_rom),max(f)]) + 0.03;
        min_x_store(i) = min_x;
        max_x_store(i) = max_x;
        max_y_store(i) = max_y;
%     else
%         min_x = min_x_store(i);
%         max_x = max_x_store(i);
%         max_y = max_y_store(i);
%     end
    xlim([min_x max_x]);    
    ylim([0 max_y]);
%     legend('Truth','ture coeffs','ROM','location','northeast');
%     legend('Truth','Linear model','location','northwest');
%     legend('boxoff');

    set(gca,'fontsize',12,'TickLabelInterpreter','latex');
%     set(gca,'YScale','Log');
    xlabel('TMX','fontsize',15,'interpreter','latex');
    ylabel('PDF','fontsize',15);
    lat_tmp = round(city_coords.lat(i)*10)/10;
    lon_tmp = round(city_coords.lon(i)*10)/10;
    title_name = strcat(city_coords.city(i),{' '},num2str(lat_tmp),'N,',...
        num2str(lon_tmp),'E',{', '},num2str(year,'%04d'),...
        '-',num2str(year+year_skip-1,'%04d'));
    title(title_name);
    set(gcf,'Units','inches',...
        'Position',[1 1 4.2 3],'PaperPosition',[0 0 4.2 3]);
    if(flag_print == 1)
        outname = strcat('figs/0207_localPDF_',num2str(i,'%02d'),'_',num2str(year,'%04d'));
        print(outname,'-djpeg','-r200');
        close;
    end
end
end
%%
elseif(choose_post == 2)
load coastlines;
% loop over location
for i = 7 %1:40

    data_true = [];
    data_est = [];
    data_rom = [];
    data_tm = [];

    for year = 1979:year_skip:1998
        filename = strcat('stats_4year/localPDF_4year_',num2str(year,'%d'),'.mat');
        load(filename);
        tmp = Tmx_local(:,i) - 273.15;
        data_true = [data_true; tmp];
        filename = strcat('stats_4year/localPDF_est_4year_',num2str(year,'%d'),'.mat');
        load(filename);
        tmp = Tmx_local(:,i) - 273.15;
        data_est = [data_est; tmp];
        for kk=1:10
            filename = strcat('stats_4year/localPDF_rom_4year_',num2str(year,'%d'),...,
                '_',num2str(kk,'%03d'),'.mat');
            load(filename);
            in_buffer(:,kk) = Tmx_local(:,i);
        end
        tmp = in_buffer - 273.15;
        data_rom = [data_rom; tmp];
    end

    [f,xi] = ksdensity(data_true);
    % bootstrapping true PDF
    epsilon = 1e-8;
    nboot = 10;
    pdf_boot = zeros(nboot,numel(xi));
    for k=1:nboot
        sample = datasample(data_true,numel(data_true));
        pdf_boot(k,:) = ksdensity(sample,xi);
    end
    pdf_error = std(pdf_boot);
    f_low = f - pdf_error;
    f_upp = f + pdf_error;
    f_low(f_low < 0) = epsilon;
    %     logpdf_error = std(log10(pdf_boot));
    %     f_low = f./(10.^logpdf_error);
    %     f_upp = f.*(10.^logpdf_error);

    [f_est,xi_est] = ksdensity(data_est);
    [f_rom,xi_rom] = ksdensity(data_rom(:));

    % PDF of ROM 
    f_rom_mrlz = zeros(numel(xi_rom),10);
    for l=1:10
        f_rom_mrlz(:,l) = ksdensity(data_rom(:,l),xi_rom);
    end
    f_rom_error = std(f_rom_mrlz,[],2)';
    f_rom_low = f_rom - f_rom_error;
    f_rom_upp = f_rom + f_rom_error;
    f_rom_low(f_rom_low < 0) = epsilon;
    %     f_tm_low = zeros(size(f_tm));
    %     f_tm_upp = zeros(size(f_tm));
    %     for j = 1:numel(xi_tm)
    %         f_tm_low(j) = min(f_tm_mrlz(j,:));
    %         f_tm_upp(j) = max(f_tm_mrlz(j,:));
    %     end

    % plot
    figure;
    fill([xi,fliplr(xi)],[f_low,fliplr(f_upp)],...
        [1 0.7 0.7],'EdgeColor','none');
    hold on;
    fill([xi_rom,fliplr(xi_rom)],[f_rom_low,fliplr(f_rom_upp)],...
        [0.7 0.7 0.7],'EdgeColor','none');
    plot(xi,f,'r--','linewidth',2);
    hold on;

    plot(xi_rom,f_rom,'k','linewidth',2);
    min_x = min([min(xi_rom),min(xi)]);
    max_x = max([max(xi_rom),max(xi)]);
    xlim([min_x, max_x]);
    max_y = max([max(f_rom),max(f)]);
    %     ylim([1e-8 1]);

    set(gca,'fontsize',12,'TickLabelInterpreter','latex');
%     set(gca,'YScale','Log');
    xlabel('TMX','fontsize',15,'interpreter','latex');
    ylabel('PDF','fontsize',15);
    lat_tmp = round(city_coords.lat(i)*10)/10;
    lon_tmp = round(city_coords.lon(i)*10)/10;
    title_name = strcat(city_coords.city(i),{' '},num2str(lat_tmp),'N,',...
        num2str(lon_tmp),'E');
    title(title_name);
    set(gcf,'Units','inches',...
        'Position',[1 1 4.2 3]);
    %     outname = strcat('0822_US_20year_localPDF_linear_',num2str(i,'%02d'));
    %     print(outname,'-djpeg','-r300');
    %     close;
end
elseif(choose_post == 3)
load coastlines;
load('UScont/UScont_mode_maximum.mat');
% loop over location
for i = 1:10 %1:40

    data_true = [];
    data_est = [];
    data_rom = [];
    data_tm = [];

    for year = 1979:year_skip:1998
        filename = strcat('UScont/localPDF_modemax_US_4year_',num2str(year,'%d'),'.mat');
        load(filename);
        tmp = Tmx_local(:,i) - 273.15;
        data_true = [data_true; tmp];
        filename = strcat('UScont/localPDF_modemax_US_est_4year_',num2str(year,'%d'),'.mat');
        load(filename);
        tmp = Tmx_local(:,i) - 273.15;
        data_est = [data_est; tmp];
        filename = strcat('UScont/localPDF_modemax_US_rom_4year_',num2str(year,'%d'),'.mat');
        load(filename);
        tmp = Tmx_local(:,i,:) - 273.15;
        tmp = tmp(:);
        data_rom = [data_rom; tmp];
        filename = strcat('UScont/localPDF_modemax_US_rom_atm_4year_',num2str(year,'%d'),'.mat');
        load(filename);
        tmp = Tmx_local(:,i,:) - 273.15;
        tmp = squeeze(tmp);
        %     tmp = tmp(:);
        data_tm = [data_tm; tmp];
    end

    [f,xi] = ksdensity(data_true);
    % bootstrapping true PDF
    epsilon = 1e-8;
    nboot = 10;
    pdf_boot = zeros(nboot,numel(xi));
    for k=1:nboot
        sample = datasample(data_true,numel(data_true));
        pdf_boot(k,:) = ksdensity(sample,xi);
    end
    pdf_error = std(pdf_boot);
    f_low = f - pdf_error;
    f_upp = f + pdf_error;
    f_low(f_low < 0) = epsilon;
    %     logpdf_error = std(log10(pdf_boot));
    %     f_low = f./(10.^logpdf_error);
    %     f_upp = f.*(10.^logpdf_error);

    [f_est,xi_est] = ksdensity(data_est);
    [f_rom,xi_rom] = ksdensity(data_rom);

    % PDF of ROM with TM
    [f_tm,xi_tm] = ksdensity(data_tm(:));
    f_rom_mrlz = zeros(numel(xi_tm),10);
    for l=1:10
        f_rom_mrlz(:,l) = ksdensity(data_tm(:,l),xi_tm);
    end
    f_rom_error = std(f_rom_mrlz,[],2)';
    f_rom_low = f_tm - f_rom_error;
    f_rom_upp = f_tm + f_rom_error;
    f_rom_low(f_rom_low < 0) = epsilon;
    %     f_tm_low = zeros(size(f_tm));
    %     f_tm_upp = zeros(size(f_tm));
    %     for j = 1:numel(xi_tm)
    %         f_tm_low(j) = min(f_tm_mrlz(j,:));
    %         f_tm_upp(j) = max(f_tm_mrlz(j,:));
    %     end

    % plot
    figure;
    fill([xi,fliplr(xi)],[f_low,fliplr(f_upp)],...
        [1 0.7 0.7],'EdgeColor','none');
    hold on;
    fill([xi_tm,fliplr(xi_tm)],[f_rom_low,fliplr(f_rom_upp)],...
        [0.7 0.7 1],'EdgeColor','none');
    plot(xi,f,'r--','linewidth',2);
    hold on;

    plot(xi_rom,f_rom,'k','linewidth',2);

    plot(xi_tm,f_tm,'b','linewidth',2);
    min_x = min([min(xi_rom),min(xi_tm),min(xi)]);
    max_x = max([max(xi_rom),max(xi_tm),max(xi)]);
    xlim([min_x, max_x]);
    max_y = max([max(f_rom),max(f)]);
    %     ylim([1e-8 1]);

    set(gca,'fontsize',12,'TickLabelInterpreter','latex');
%     set(gca,'YScale','Log');
    xlabel('TMX','fontsize',15,'interpreter','latex');
    ylabel('PDF','fontsize',15);
    title_name = strcat(num2str(lat_max(i)),'N,',num2str(lon_max(i)-360),'W');
    title(title_name);
    set(gcf,'Units','inches',...
        'Position',[1 1 4.2 3]);
        outname = strcat('0828_US_20year_localPDF_modemax_linear_',num2str(i,'%02d'));
        print(outname,'-djpeg','-r300');
        close;
end
end