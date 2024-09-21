clear;
close all;

choose_post = 6; % 1 for TMX averaged over 10 years
                 % 2 for std of TMX
                 % 2.5 for std, northern hemisphere
                 % 6 for quantiles, 10 years
                 % 6.5 for quantiles, northern hemisphere
                 % 7 for quantiles, subtract 79-82 mean
                 % 8 for quantiles, 20 years

year_skip = 10;  % number of years sampled altogether

n_ens = 25;

lat_skip = 1;
lon_skip = 1;
%% load grid and get the coarse grid
folder = 'download/';
variable = 'tasmax';

load('pca_historical/Tmx_clim_historical.mat','lon','lat','Tmx_mean');
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
n_ens_stats = 50;
for year = 1850:year_skip:1890
    for l=1:n_ens_stats
        filename = strcat('stats_10year/stats_10year_',scenario,num2str(year,'%d'),'_',num2str(l,'%03d'),'.mat');
        load(filename);
        Tmx_mean_ref = Tmx_mean_ref + TMX_avg;
        Tmx_std_ref = Tmx_std_ref + TMX_std.^2 + TMX_avg.^2;
        filename = strcat('stats_10year/qntl_10year_',scenario,num2str(year,'%d'),'_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_ref = q_ref + qntl_975;
    end
    n = n + 1;
end
Tmx_mean_ref = Tmx_mean_ref / n / n_ens_stats;
Tmx_std_ref = Tmx_std_ref / n / n_ens_stats;
Tmx_std_ref = Tmx_std_ref - Tmx_mean_ref.^2;
Tmx_std_ref = sqrt(Tmx_std_ref);
q_ref = q_ref / n/ n_ens_stats;
file_output = strcat('stats_10year/stats_historical_1850_1899_',num2str(n_ens_stats,'%d'),'.mat');
save(file_output,'Tmx_mean_ref','Tmx_std_ref','q_ref','n_ens_stats');

end
%%
load('stats_10year/stats_historical_1850_1899_50.mat','Tmx_mean_ref','Tmx_std_ref','q_ref');
if(choose_post == 1)
    load coastlines;
    n_rlz = 50;
    levels_Tmx = linspace(-2,10,20);
    levels = linspace(-2,2,10);
%     scenario = '';
%     for year = 2000%:year_skip:2000
    scenario = 'ssp126_';
    for year = 2090%:year_skip:2090

        Tmx_mean_true = zeros(n_lonc,n_latc);
        for l=1:n_ens
            filename = strcat('stats_10year/ref/stats_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_true = Tmx_mean_true + TMX_avg;
        end
        Tmx_mean_true = Tmx_mean_true / n_ens;
        

        Tmx_mean_rom = zeros(size(Tmx_mean_true));
        Tmx_mean_rom_mrlz = zeros(n_lonc,n_latc,n_rlz);
        for l = 1:n_rlz
            filename = strcat('stats_10year/rom_ens_50/stats_rom_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_rom = Tmx_mean_rom + TMX_avg;
            Tmx_mean_rom_mrlz(:,:,l) = TMX_avg;
        end
        Tmx_mean_rom = Tmx_mean_rom / n_rlz;

        Tmx_mean_ml = zeros(size(Tmx_mean_true));
        Tmx_mean_ml_mrlz = zeros(n_lonc,n_latc,n_rlz);
        for l = 1:n_rlz
            filename = strcat('stats_10year/rom_cov_ens_50/stats_rom_cov_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_ml = Tmx_mean_ml + TMX_avg;
            Tmx_mean_ml_mrlz(:,:,l) = TMX_avg;
        end
        Tmx_mean_ml = Tmx_mean_ml / n_rlz;

        % compute anmolay from the reference period
        Tmx_mean_true = Tmx_mean_true - Tmx_mean_ref;
        Tmx_mean_rom = Tmx_mean_rom - Tmx_mean_ref;
        Tmx_mean_ml = Tmx_mean_ml - Tmx_mean_ref;

        idx = isnan(Tmx_mean_rom);
        Tmx_mean_true(idx) = NaN;

        lon_plot = [lonc(n_lonc/2+1:n_lonc) - 360;lonc(1:n_lonc/2)];
        mu_plot = [Tmx_mean_true(n_lonc/2+1:n_lonc,:); Tmx_mean_true(1:n_lonc/2,:)];

        figure;
        set(gcf,'Units','inches',...
            'Position',[1 1 20 6],'PaperPosition',[0 0 21 7]);

        h1 = subplot('Position',[0.05 0.05 0.27 0.9]);
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels_Tmx,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.0500    0.1531    0.27    0.7969];
        colormap(h1,flipud(hot));
        caxis([min(levels_Tmx) max(levels_Tmx)]);
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
% % 
        h3 = subplot('Position',[0.35 0.05 0.27 0.9]);
%         mu_plot = [Tmx_mean_rom(n_lonc/2+1:n_lonc,:); Tmx_mean_rom(1:n_lonc/2,:)];
%         contourf(lon_plot,latc,mu_plot' - 273.15,20,'linestyle','none');
        error_rom = Tmx_mean_rom(:,:) - Tmx_mean_true(:,:);
        mu_plot = [error_rom(n_lonc/2+1:n_lonc,:); error_rom(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.3500    0.1531    0.27    0.7969];
        colormap(h3,redblue);
        caxis([-2 2]);
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
        error_ml = Tmx_mean_ml(:,:) - Tmx_mean_true(:,:);
        mu_plot = [error_ml(n_lonc/2+1:n_lonc,:); error_ml(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.6500    0.1531    0.27    0.7969];
        colormap(h4,redblue);
        caxis([-2 2]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        mlabel('off');
        set(gca,'TickLabelInterpreter', 'latex');
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);

        
%         sgtitle(strcat(num2str(year,'%d'),'-',num2str(year+year_skip-1,'%d')));

        outname = strcat('figs_manu/MPI_ens50_',scenario,'mean_',num2str(year,'%d'));
        print(outname,'-djpeg','-r300');
%         close;
    end
elseif(choose_post == 2)
    load coastlines;
    levels = linspace(-2,2,10);
    n_rlz = 50;
%     scenario = '';   
%     for year = 2000 %:year_skip:2000
    scenario = 'ssp585_';
    for year = 2090%:year_skip:2090

        Tmx_std_true = zeros(n_lonc,n_latc);
        Tmx_mean_true = zeros(n_lonc,n_latc);
        for l=1:n_ens

            filename = strcat('stats_10year/ref/stats_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_true = Tmx_mean_true + TMX_avg;
            Tmx_std_true = Tmx_std_true + TMX_std.^2 + TMX_avg.^2;
        end
        Tmx_mean_true = Tmx_mean_true / n_ens;
        Tmx_std_true = Tmx_std_true / n_ens;
        Tmx_std_true = Tmx_std_true - Tmx_mean_true.^2;
        Tmx_std_true = sqrt(Tmx_std_true);

        Tmx_std_rom = zeros(size(Tmx_std_true));
        Tmx_std_rom_mrlz = zeros(n_lonc,n_latc,10);
        Tmx_mean_rom = zeros(size(Tmx_std_true));
        for l = 1:n_rlz
            filename = strcat(['stats_10year/rom_ens_50/' ...
                'stats_rom_10year_'],scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_rom = Tmx_mean_rom + TMX_avg;
            Tmx_std_rom = Tmx_std_rom + TMX_std.^2 + TMX_avg.^2;
            Tmx_std_rom_mrlz(:,:,l) = TMX_std;
        end
        Tmx_mean_rom = Tmx_mean_rom / n_rlz;
        Tmx_std_rom = Tmx_std_rom / n_rlz;
        Tmx_std_rom = Tmx_std_rom - Tmx_mean_rom.^2;
        Tmx_std_rom = sqrt(Tmx_std_rom);

        % load COV ROM results
        Tmx_std_ml = zeros(size(Tmx_std_true));
        Tmx_std_ml_mrlz = zeros(n_lonc,n_latc,10);
        Tmx_mean_ml = zeros(size(Tmx_std_true));
        for l = 1:n_rlz
            filename = strcat(['stats_10year/rom_cov_ens_50/' ...
                'stats_rom_cov_10year_'],scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            Tmx_mean_ml = Tmx_mean_ml + TMX_avg;
            Tmx_std_ml = Tmx_std_ml + TMX_std.^2 + TMX_avg.^2;
            Tmx_std_ml_mrlz(:,:,l) = TMX_std;
        end
        Tmx_mean_ml = Tmx_mean_ml / n_rlz;
        Tmx_std_ml = Tmx_std_ml / n_rlz;
        Tmx_std_ml = Tmx_std_ml - Tmx_mean_ml.^2;
        Tmx_std_ml = sqrt(Tmx_std_ml);

%         Tmx_std_true = Tmx_std_true - Tmx_std_ref;
%         Tmx_std_rom_noice = Tmx_std_rom_noice - Tmx_std_ref;
        idx = isnan(Tmx_std_rom);
        Tmx_std_true(idx) = NaN;
        lon_plot = [lonc(n_lonc/2+1:n_lonc) - 360;lonc(1:n_lonc/2)];
        mu_plot = [Tmx_std_true(n_lonc/2+1:n_lonc,:); Tmx_std_true(1:n_lonc/2,:)];

        figure;
        set(gcf,'Units','inches',...
            'Position',[1 1 20 6],'PaperPosition',[0 0 21 7]);

        h1 = subplot('Position',[0.05 0.05 0.27 0.9]);
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',linspace(0,8,20),'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.0500    0.1531    0.27    0.7969];
        colormap(h1,jet);
        caxis([0 8]);
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
        error_rom = Tmx_std_rom - Tmx_std_true;
        mu_plot = [error_rom(n_lonc/2+1:n_lonc,:); error_rom(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        mlabel('off');
        cb = colorbar('location','southoutside');
        ax_pos = [0.3500    0.1531    0.27    0.7969];
        colormap(h3,redblue);
        caxis([-2 2]);
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
        error_ml = Tmx_std_ml - Tmx_std_true;
        mu_plot = [error_ml(n_lonc/2+1:n_lonc,:); error_ml(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        mlabel('off');
        cb = colorbar('location','southoutside');
        ax_pos = [0.6500    0.1531    0.27    0.7969];
        colormap(h4,redblue);
        caxis([-2 2]);
        set(gca,'fontsize',12);
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        xlabel('longitude','fontsize',15);
        ylabel('latitude','fontsize',15);
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.02]);
        set(ax,'Position',ax_pos);
        set(gca,'TickLabelInterpreter', 'latex');

        sgtitle(strcat(num2str(year,'%d'),'-',num2str(year+year_skip-1,'%d')));
        

        outname = strcat('figs_manu/MPI_ens50_',scenario,'std_',num2str(year,'%d'));
        print(outname,'-djpeg','-r300');
%         close;
    end
elseif(choose_post == 6)
    load coastlines;
    levels = linspace(-3,3,10);
    levels_Tmx = linspace(-2,10,20);
    n_rlz = 50;
%     scenario = '';
%     for year = 2000%:year_skip:2000
    scenario = 'ssp126_';
    for year = 2090%:year_skip:2090

        q_true = zeros(n_lonc,n_latc);
        for l=1:n_ens

            filename = strcat('stats_10year/ref/qntl_10year_',scenario,num2str(year,'%d'),...
                '_',num2str(l,'%03d'),'.mat');
            load(filename);
            q_true = q_true + qntl_975;
        end
        q_true = q_true / n_ens;


        q_rom_mrlz = zeros(n_lonc,n_latc,n_rlz);
        for l=1:n_rlz
        filename = strcat('stats_10year/rom_ens_50/qntl_rom_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_rom_mrlz(:,:,l) = qntl_975;
        end
        q_rom = mean(q_rom_mrlz,3);
        
        q_ml_mrlz = zeros(n_lonc,n_latc,n_rlz);
        for l=1:n_rlz
        filename = strcat('stats_10year/rom_cov_ens_50/qntl_rom_cov_10year_',scenario,num2str(year,'%d'),...
            '_',num2str(l,'%03d'),'.mat');
        load(filename);
        q_ml_mrlz(:,:,l) = qntl_975;
        end
        q_ml = mean(q_ml_mrlz,3);

        q_true = q_true - q_ref;
        q_rom = q_rom - q_ref;
        q_ml = q_ml - q_ref;
        idx = isnan(q_rom);
        q_true(idx) = NaN;

        lon_plot = [lonc(n_lonc/2+1:n_lonc) - 360;lonc(1:n_lonc/2)];

        set(gcf,'Units','inches',...
            'Position',[1 1 20 6],'PaperPosition',[0 0 21 7]);

        h1 = subplot('Position',[0.05 0.05 0.27 0.9]);
        mu_plot = [q_true(n_lonc/2+1:n_lonc,:); q_true(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels_Tmx,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.0500    0.1531    0.27    0.7969];
        mlabel('off');
        colormap(h1,flipud(hot));
        caxis([min(levels_Tmx) max(levels_Tmx)]);
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
        mu_plot = [q_rom(n_lonc/2+1:n_lonc,:); q_rom(1:n_lonc/2,:)] - ...
                  [q_true(n_lonc/2+1:n_lonc,:); q_true(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.3500    0.1531    0.27    0.7969];
        mlabel('off');
        colormap(h3,redblue);
        caxis([-3 3]);
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
        mu_plot = [q_ml(n_lonc/2+1:n_lonc,:); q_ml(1:n_lonc/2,:)] - ...
                  [q_true(n_lonc/2+1:n_lonc,:); q_true(1:n_lonc/2,:)];
        ax = worldmap('world');
        contourfm(latc,lon_plot,mu_plot',levels,'linestyle','none');
        hold on;
        geoshow(coastlat,coastlon,'color','k');
        cb = colorbar('location','southoutside');
        ax_pos = [0.6500    0.1531    0.27    0.7969];
        mlabel('off');
        colormap(h4,redblue);
        caxis([-3 3]);
        set(gca,'fontsize',12);
        set(gca,'TickLabelInterpreter', 'latex');
        setm(ax,'FFaceColor', [0.7 0.7 0.7]);
        cb_pos = cb.Position;
        set(cb,'Position',[cb_pos(1)+0.01,cb_pos(2)+0.01,cb_pos(3)-0.02,cb_pos(4)-0.01]);
        set(ax,'Position',ax_pos);   

        
        sgtitle(strcat(num2str(year,'%d'),'-',num2str(year+year_skip-1,'%d')));
        

        outname = strcat('figs_manu/MPI_ens50_',scenario,'qntl_',num2str(year,'%d'));
        print(outname,'-djpeg','-r300');
%         close;
    end
end