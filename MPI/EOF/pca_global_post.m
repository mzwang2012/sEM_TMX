clear;
close all;

load('pca_historical/Tmx_clim_historical.mat','lon','lat');

n_lat = length(lat);
n_lon = length(lon);

%%
load('pca_historical/modes_Tmx_365skip5_global_historical_ens.mat','w_hat');
mode = w_hat;
clear w_hat;

%%
% load('pca_historical/modes_Tmx_365skip5_global_historical_ens.mat','lambda');
% 
% ratio = cumsum(lambda);
% % semilogy(lambda(1:end-365)/lambda(1),'k','linewidth',2);
% plot(ratio / ratio(end),'k','linewidth',2);
% xlim([0 5000]);
% ylim([0.8 1]);
% xlabel('mode','fontsize',15);
% % ylabel('$\lambda / \lambda_1$','interpreter','latex','fontsize',15);
% ylabel('ratio of energy','fontsize',15);
% set(gca,'TickLabelInterpreter', 'latex','fontsize',12);
% set(gcf,'Units','inches',...
%     'Position',[1 1 4.5 2.5]);
%% plot
for i_mode = 3

lon_plot = [lon(n_lon/2+1:n_lon) - 360;lon(1:n_lon/2)];
true_plot = [mode(n_lon/2+1:n_lon,:,i_mode); mode(1:n_lon/2,:,i_mode)];

load coastlines;

figure;
% contourf(lon_plot,lat,true_plot',20,'linestyle','none');
% hold on;
% plot(coastlon,coastlat,'k');
worldmap('world');
contourfm(lat,lon_plot,true_plot',20,'linestyle','none');
hold on;
geoshow(coastlat,coastlon,'color','k');
colorbar;
colormap(redblue);
mlabel('off');
% caxis([-1e-2 1e-2]);
% caxis([0.8*min(true_plot(:)),0.8*max(true_plot(:))]);
set(gca,'fontsize',12);
daspect([1 1 1]);
% title_name = strcat('rank',{' '},num2str(kk,'%d'),{', '},...
%     'mode',{' '},num2str(i_mode,'%d'));
%         title(title_name);
% xlabel('longitude','fontsize',15);
% ylabel('latitude','fontsize',15);
% xlim([-180 180]);
% ylim([-90 90]);
set(gca,'TickLabelInterpreter', 'latex');
set(gcf,'Units','inches',...
    'Position',[1 1 8 3.5],'PaperPosition',[0 0 9 3.5]);
% outname = strcat('figs/1204_pca_hist_mode_',num2str(i_mode,'%3d'));
% print(outname,'-djpeg','-r200');
% close;
end