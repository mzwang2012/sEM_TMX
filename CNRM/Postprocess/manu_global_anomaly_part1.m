clear;
close all;

choose_post = 1; % 1 for TMX averaged over 10 years
                 % 2 for std of TMX
                 % 2.5 for std, northern hemisphere
                 % 6 for quantiles, 10 years
                 % 6.5 for quantiles, northern hemisphere
                 % 7 for quantiles, subtract 79-82 mean
                 % 8 for quantiles, 20 years

year_skip = 10;  % number of years sampled altogether

lat_skip = 2;
lon_skip = 2;
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

%% compute and store historical average
if(1==0)
scenario = '';

Tmx_mean_ref = zeros(n_lonc,n_latc);
Tmx_std_ref = zeros(n_lonc,n_latc);
q_ref = zeros(n_lonc,n_latc);
n = 0;
for year = 1850:year_skip:1890
    filename = strcat('stats_10year/stats_10year_',scenario,num2str(year,'%d'),'.mat');
    load(filename);
    Tmx_mean_ref = Tmx_mean_ref + TMX_avg;
    Tmx_std_ref = Tmx_std_ref + TMX_std.^2 + TMX_avg.^2;
    filename = strcat('stats_10year/qntl_10year_',scenario,num2str(year,'%d'),'.mat');
    load(filename);
    q_ref = q_ref + qntl_975;
    n = n + 1;
end
Tmx_mean_ref = Tmx_mean_ref / n;
Tmx_std_ref = Tmx_std_ref / n;
Tmx_std_ref = Tmx_std_ref - Tmx_mean_ref.^2;
Tmx_std_ref = sqrt(Tmx_std_ref);
q_ref = q_ref / n;
save('stats_noice/stats_historical_1850_1899.mat','Tmx_mean_ref','Tmx_std_ref','q_ref');

end
%%
load('stats_noice/stats_historical_1850_1899.mat','Tmx_mean_ref','Tmx_std_ref','q_ref');
% mean Tmx(lon,lat,year) 
if(choose_post == 1)
    load coastlines;
    n_rlz = 10;
    levels = linspace(-2,10,20);
%     scenario = '';
%     for year = 2000%:year_skip:2000
    scenario = 'ssp585_';
    for year = 2090%:year_skip:2090

        filename = strcat('stats_10year/stats_10year_',scenario,num2str(year,'%d'),'.mat');
        load(filename);
        Tmx_mean_true = TMX_avg;
        
%         Tmx_mean_rom_global = zeros(size(Tmx_mean_true));
%         for l = 1:n_rlz
%             filename = strcat('stats_10year/rom_glbT_lin_mean_var/stats_rom_10year_',scenario,num2str(year,'%d'),...
%                 '_',num2str(l,'%03d'),'.mat');
%             load(filename);
%             Tmx_mean_rom_global = Tmx_mean_rom_global + TMX_avg;
%         end
%         Tmx_mean_rom_global = Tmx_mean_rom_global / n_rlz;

        Tmx_mean_rom_noice = zeros(size(Tmx_mean_true));
        Tmx_mean_rom_noice_mrlz = zeros([size(Tmx_mean_true),n_rlz]);
        for l = 1:n_rlz
            filename = strcat('stats_10year/stats_rom_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_rom_noice = Tmx_mean_rom_noice + TMX_avg;
            Tmx_mean_rom_noice_mrlz(:,:,l) = TMX_avg;
        end
        Tmx_mean_rom_noice = Tmx_mean_rom_noice / n_rlz;

        Tmx_mean_ml_noice = zeros(size(Tmx_mean_true));
        for l = 1:n_rlz
            filename = strcat('stats_10year/stats_rom_ML_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_ml_noice = Tmx_mean_ml_noice + TMX_avg;
        end
        Tmx_mean_ml_noice = Tmx_mean_ml_noice / n_rlz;

        % compute anmolay from the reference period
        Tmx_mean_true = Tmx_mean_true - Tmx_mean_ref;
%         Tmx_mean_rom_global = Tmx_mean_rom_global - Tmx_mean_ref;
        Tmx_mean_rom_noice = Tmx_mean_rom_noice - Tmx_mean_ref;
        Tmx_mean_rom_noice_mrlz = Tmx_mean_rom_noice_mrlz - Tmx_mean_ref;
        Tmx_mean_ml_noice = Tmx_mean_ml_noice - Tmx_mean_ref;

        idx = isnan(Tmx_mean_rom_noice);
        Tmx_mean_true(idx) = NaN;
%         Tmx_mean_rom_global(idx) = NaN;

        lon_plot = [lonc(n_lonc/2+1:n_lonc) - 360;lonc(1:n_lonc/2)];
        mu_plot = [Tmx_mean_true(n_lonc/2+1:n_lonc,:); Tmx_mean_true(1:n_lonc/2,:)];

        figure;
        set(gcf,'Units','inches',...
            'Position',[1 1 20 6],'PaperPosition',[0 0 21 7]);

        h1 = subplot('Position',[0.05 0.05 0.27 0.9]);
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.0500    0.1531    0.27    0.7969];
        colormap(h1,flipud(hot));
        caxis([min(levels) max(levels)]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        mlabel('off');
        set(gca,'TickLabelInterpreter', 'latex');
        % Adjust size of colorbar and contour plot
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);
%         title('True field');
        
% 
        h3 = subplot('Position',[0.35 0.05 0.27 0.9]);
%         mu_plot = [Tmx_mean_rom(n_lonc/2+1:n_lonc,:); Tmx_mean_rom(1:n_lonc/2,:)];
%         contourf(lon_plot,latc,mu_plot' - 273.15,20,'linestyle','none');
        error_rom = Tmx_mean_rom_noice_mrlz(:,:,1); %- Tmx_mean_true(:,:);
        mu_plot = [error_rom(n_lonc/2+1:n_lonc,:); error_rom(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.3500    0.1531    0.27    0.7969];
        colormap(h3,flipud(hot));
        caxis([min(levels) max(levels)]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        mlabel('off');
        set(gca,'TickLabelInterpreter', 'latex');
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);
%         title('ROM without ice (error)');

        h4 = subplot('Position',[0.65 0.05 0.28 0.9]);
%         mu_plot = [Tmx_mean_rom(n_lonc/2+1:n_lonc,:); Tmx_mean_rom(1:n_lonc/2,:)];
%         contourf(lon_plot,latc,mu_plot' - 273.15,20,'linestyle','none');
        error_ml = Tmx_mean_rom_noice_mrlz(:,:,2); % - Tmx_mean_true(:,:);
        mu_plot = [error_ml(n_lonc/2+1:n_lonc,:); error_ml(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.6500    0.1531    0.27    0.7969];
        colormap(h4,flipud(hot));
        caxis([min(levels) max(levels)]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        mlabel('off');
        set(gca,'TickLabelInterpreter', 'latex');
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);

        
%         sgtitle(strcat(num2str(year,'%d'),'-',num2str(year+year_skip-1,'%d')));

        outname = strcat('figs_manu/global_',scenario,'mean_',num2str(year,'%d'),'_part1');
        print(outname,'-djpeg','-r300');
%         close;
    end
elseif(choose_post == 2)
    load coastlines;
    levels = linspace(0,8,20);
    n_rlz = 10;
    scenario = '';   
    for year = 2000 %:year_skip:2000
%     scenario = 'ssp126_';
%     for year = 2090%:year_skip:2090

        filename = strcat('stats_10year/stats_10year_',scenario,num2str(year,'%d'),'.mat');
        load(filename);
        Tmx_std_true = TMX_std;
%         filename = strcat('stats_10year/stats_est_10year_',scenario,num2str(year,'%d'),'.mat');
%         filename = strcat('stats_10year/rom_glbT_true_mean_var/',...
%             'stats_rom_10year_',scenario,num2str(year,'%d'),'_001.mat');
%         load(filename);
%         Tmx_std_est = TMX_std;

%         Tmx_std_rom_global = zeros(size(Tmx_std_true));
% %         Tmx_std_rom_mrlz = zeros(n_lonc,n_latc,10);
%         Tmx_mean_rom_global = zeros(size(Tmx_std_true));
%         for l = 1:n_rlz
%             filename = strcat(['stats_10year/rom_glbT_lin_mean_var/' ...
%                 'stats_rom_10year_'],scenario,num2str(year,'%d'),...
%                 '_',num2str(l,'%03d'),'.mat');
%             load(filename);
%             Tmx_mean_rom_global = Tmx_mean_rom_global + TMX_avg;
%             Tmx_std_rom_global = Tmx_std_rom_global + TMX_std.^2 + TMX_avg.^2;
% %             Tmx_std_rom_mrlz(:,:,l) = TMX_std;
%         end
%         Tmx_mean_rom_global = Tmx_mean_rom_global / n_rlz;
%         Tmx_std_rom_global = Tmx_std_rom_global / n_rlz;
%         Tmx_std_rom_global = Tmx_std_rom_global - Tmx_mean_rom_global.^2;
%         Tmx_std_rom_global = sqrt(Tmx_std_rom_global);

        % load no ice ROM results
        Tmx_std_rom_noice = zeros(size(Tmx_std_true));
        Tmx_mean_rom_noice = zeros(size(Tmx_std_true));
        Tmx_std_rom_noice_mrlz = zeros(n_lonc,n_latc,10);
        for l = 1:n_rlz
            filename = strcat(['stats_10year/' ...
                'stats_rom_10year_'],scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_rom_noice = Tmx_mean_rom_noice + TMX_avg;
            Tmx_std_rom_noice = Tmx_std_rom_noice + TMX_std.^2 + TMX_avg.^2;
            Tmx_std_rom_noice_mrlz(:,:,l) = TMX_std;
        end
        Tmx_mean_rom_noice = Tmx_mean_rom_noice / n_rlz;
        Tmx_std_rom_noice = Tmx_std_rom_noice / n_rlz;
        Tmx_std_rom_noice = Tmx_std_rom_noice - Tmx_mean_rom_noice.^2;
        Tmx_std_rom_noice = sqrt(Tmx_std_rom_noice);

        % load no ice ML results
        Tmx_std_ml_noice = zeros(size(Tmx_std_true));
        Tmx_mean_ml_noice = zeros(size(Tmx_std_true));
        for l = 1:n_rlz
            filename = strcat(['stats_10year/' ...
                'stats_rom_ML_10year_'],scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_ml_noice = Tmx_mean_ml_noice + TMX_avg;
            Tmx_std_ml_noice = Tmx_std_ml_noice + TMX_std.^2 + TMX_avg.^2;
        end
        Tmx_mean_ml_noice = Tmx_mean_ml_noice / n_rlz;
        Tmx_std_ml_noice = Tmx_std_ml_noice / n_rlz;
        Tmx_std_ml_noice = Tmx_std_ml_noice - Tmx_mean_ml_noice.^2;
        Tmx_std_ml_noice = sqrt(Tmx_std_ml_noice); 

%         Tmx_std_true = Tmx_std_true - Tmx_std_ref;
%         Tmx_std_rom_noice = Tmx_std_rom_noice - Tmx_std_ref;
        idx = isnan(Tmx_std_rom_noice);
        Tmx_std_true(idx) = NaN;
        lon_plot = [lonc(n_lonc/2+1:n_lonc) - 360;lonc(1:n_lonc/2)];
        mu_plot = [Tmx_std_true(n_lonc/2+1:n_lonc,:); Tmx_std_true(1:n_lonc/2,:)];

        figure;
        set(gcf,'Units','inches',...
            'Position',[1 1 20 6],'PaperPosition',[0 0 21 7]);

        h1 = subplot('Position',[0.05 0.05 0.27 0.9]);
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.0500    0.1531    0.27    0.7969];
        colormap(h1,jet);
        caxis([min(levels),max(levels)]);
        mlabel('off');
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        set(gca,'fontsize',12);
        xlabel('longitude','fontsize',15);
        ylabel('latitude','fontsize',15);
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);
        set(gca,'TickLabelInterpreter', 'latex');
%         title('True field');
        
       

%         h2 = subplot('Position',[0.55 0.55 0.4 0.37]);
% %         mu_plot = [Tmx_std_rom(n_lonc/2+1:n_lonc,:); Tmx_std_rom(1:n_lonc/2,:)];
%         error_rom = Tmx_std_rom_global - Tmx_std_true;
%         mu_plot = [error_rom(n_lonc/2+1:n_lonc,:); error_rom(1:n_lonc/2,:)];
%         worldmap('world');
%         contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
%         hold on;
%         geoshow(coastlat,coastlon,'color','k');
%         mlabel('off');
%         colorbar;
%         colormap(h2,redblue);
% %         caxis([0 10]);
%         caxis([-1.5 1.5]);
%         set(gca,'fontsize',12);
%         xlabel('longitude','fontsize',15);
%         ylabel('latitude','fontsize',15);
%         title('Global ROM (error)');
% 
        h3 = subplot('Position',[0.35 0.05 0.27 0.9]);
        error_rom = Tmx_std_rom_noice_mrlz(:,:,1); %- Tmx_std_true;
        mu_plot = [error_rom(n_lonc/2+1:n_lonc,:); error_rom(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        mlabel('off');
        cb = colorbar('location','southoutside');
        ax_pos = [0.3500    0.1531    0.27    0.7969];
        colormap(h3,jet);
        caxis([min(levels),max(levels)]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        xlabel('longitude','fontsize',15);
        ylabel('latitude','fontsize',15);
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);
        set(gca,'TickLabelInterpreter', 'latex');
%         title('ROM without ice (error)');

        h4 = subplot('Position',[0.65 0.05 0.27 0.9]);
        error_ml = Tmx_std_rom_noice_mrlz(:,:,2);% - Tmx_std_true;
        mu_plot = [error_ml(n_lonc/2+1:n_lonc,:); error_ml(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        mlabel('off');
        cb = colorbar('location','southoutside');
        ax_pos = [0.6500    0.1531    0.27    0.7969];
        colormap(h4,jet);
        caxis([min(levels),max(levels)]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        xlabel('longitude','fontsize',15);
        ylabel('latitude','fontsize',15);
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);
        set(gca,'TickLabelInterpreter', 'latex');

        sgtitle(strcat(num2str(year,'%d'),'-',num2str(year+year_skip-1,'%d')));
        

%         outname = strcat('figs_manu/global_',scenario,'std_',num2str(year,'%d'),'_part1');
%         print(outname,'-djpeg','-r300');
%         close;
    end
elseif(choose_post == 6)
    load coastlines;
    levels = linspace(-2,10,20);
    n_rlz = 10;
%     scenario = '';
%     for year = 2000%:year_skip:2000
    scenario = 'ssp585_';
    for year = 2090%:year_skip:2090

        filename = strcat('stats_10year/qntl_10year_',scenario,num2str(year,'%d'),'.mat');
        load(filename);
        q_true = qntl_975;
%         filename = strcat('stats_10year/qntl_est_10year_',scenario,num2str(year,'%d'),'.mat');
%         filename = strcat('stats_10year/rom_glbT_true_mean_var/qntl_rom_10year_',...
%             scenario,num2str(year,'%d'),...
%             '_',num2str(1,'%03d'),'.mat');
%         load(filename);
%         q_est = qntl_975;

%         q_rom_mrlz = zeros(n_lonc,n_latc,n_rlz);
%         for l=1:n_rlz
%         filename = strcat('stats_10year/rom_glbT_lin_mean_var/qntl_rom_10year_',scenario,num2str(year,'%d'),...
%             '_',num2str(l,'%03d'),'.mat');
%         load(filename);
%         q_rom_mrlz(:,:,l) = qntl_975;
%         end
%         q_rom_global = mean(q_rom_mrlz,3);

        q_rom_noice_mrlz = zeros(n_lonc,n_latc,n_rlz);
        for l=1:n_rlz
        filename = strcat('stats_10year/qntl_rom_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_rom_noice_mrlz(:,:,l) = qntl_975;
        end
        q_rom_noice = mean(q_rom_noice_mrlz,3);
        
        q_ml_mrlz = zeros(n_lonc,n_latc,n_rlz);
        for l=1:n_rlz
        filename = strcat('stats_10year/qntl_rom_ML_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_ml_mrlz(:,:,l) = qntl_975;
        end
        q_ml_noice = mean(q_ml_mrlz,3);

        q_true = q_true - q_ref;
        q_rom_noice = q_rom_noice - q_ref;
        q_rom_noice_mrlz = q_rom_noice_mrlz - q_ref;
        q_ml_noice = q_ml_noice - q_ref;
        idx = isnan(q_rom_noice);
        q_true(idx) = NaN;

        lon_plot = [lonc(n_lonc/2+1:n_lonc) - 360;lonc(1:n_lonc/2)];

        set(gcf,'Units','inches',...
            'Position',[1 1 20 6],'PaperPosition',[0 0 21 7]);

        h1 = subplot('Position',[0.05 0.05 0.27 0.9]);
        mu_plot = [q_true(n_lonc/2+1:n_lonc,:); q_true(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.0500    0.1531    0.27    0.7969];
        mlabel('off');
        colormap(h1,flipud(hot));
        caxis([min(levels) max(levels)]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        set(gca,'TickLabelInterpreter', 'latex');
        % set the location of colorbar
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.01]);
        set(ax,'Position',ax_pos);
%         title('True field');
        


        % part 2
%         h2 = subplot('Position',[0.55 0.55 0.4 0.37]);
% %         mu_plot = [q_rom(n_lonc/2+1:n_lonc,:); q_rom(1:n_lonc/2,:)] - 273.15;
%         mu_plot = [q_rom_global(n_lonc/2+1:n_lonc,:); q_rom_global(1:n_lonc/2,:)] - ...
%                   [q_true(n_lonc/2+1:n_lonc,:); q_true(1:n_lonc/2,:)];
%         worldmap('world');
%         contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
%         hold on;
%         geoshow(coastlat,coastlon,'color','k');
%         mlabel('off');
%         colorbar;
%         colormap(h2,redblue);
%         caxis([-3 3]); %caxis([-40 40]);
%         set(gca,'fontsize',12);
%         xlabel('longitude','fontsize',15);
%         ylabel('latitude','fontsize',15);
%         title('Global ROM (error)');
% 
%         % part 3
        h3 = subplot('Position',[0.35 0.05 0.27 0.9]);
%         mu_plot = [q_rom_mrlz(n_lonc/2+1:n_lonc,:,2); q_rom_mrlz(1:n_lonc/2,:,2)] - 273.15;
        mu_plot = [q_rom_noice_mrlz(n_lonc/2+1:n_lonc,:,1); q_rom_noice_mrlz(1:n_lonc/2,:,1)]; %- ...
                  %[q_true(n_lonc/2+1:n_lonc,:); q_true(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.3500    0.1531    0.27    0.7969];
        mlabel('off');
        colormap(h3,flipud(hot));
        caxis([min(levels) max(levels)]);
        set(gca,'fontsize',12);
        set(gca,'TickLabelInterpreter', 'latex');
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.01]);
        set(ax,'Position',ax_pos);        
%         title('ROM without ice (error)');

        % part 4
        h4 = subplot('Position',[0.65 0.05 0.27 0.9]);
%         mu_plot = [q_rom_mrlz(n_lonc/2+1:n_lonc,:,2); q_rom_mrlz(1:n_lonc/2,:,2)] - 273.15;
        mu_plot = [q_rom_noice_mrlz(n_lonc/2+1:n_lonc,:,2); q_rom_noice_mrlz(1:n_lonc/2,:,2)]; %- ...
                  %[q_true(n_lonc/2+1:n_lonc,:); q_true(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.6500    0.1531    0.27    0.7969];
        mlabel('off');
        colormap(h4,flipud(hot));
        caxis([min(levels) max(levels)]);
        set(gca,'fontsize',12);
        set(gca,'TickLabelInterpreter', 'latex');
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.01]);
        set(ax,'Position',ax_pos);   

        
        sgtitle(strcat(num2str(year,'%d'),'-',num2str(year+year_skip-1,'%d')));
        

        outname = strcat('figs_manu/global_',scenario,'qntl_',num2str(year,'%d'),'_part1');
        print(outname,'-djpeg','-r300');
%         close;
    end
end