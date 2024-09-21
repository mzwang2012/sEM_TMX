%% plot EOFs

clear;
close all;

filename = 'pca_ERA5/Tmx_365_global_1979_1998_modes.mat';
load(filename,'lat','lon');

n_lat = length(lat);
n_lon = length(lon);

%%
load('pca_historical/modes_Tmx_365skip5_global_historical.mat','w_hat');
mode = w_hat;
clear w_hat;

%%
load('pca_historical/modes_Tmx_365skip5_global_historical.mat','lambda');

ratio = cumsum(lambda);
% semilogy(lambda(1:end-365)/lambda(1),'k','linewidth',2);
plot(ratio / ratio(end),'k','linewidth',2);
hold on;
plot([0 5000],[0.95,0.95],'k--','linewidth',1.2);
ylim([0.8 1]);
xlabel('number of EOFs','fontsize',15,'fontname','Times');
ylabel('cumulative var ratio','fontsize',15,'fontname','Times');
xlim([0 5000]);
set(gca,'TickLabelInterpreter', 'latex','fontsize',12);
set(gcf,'Units','inches',...
    'Position',[1 1 4.5 2.5]);
% ylabel('$\lambda / \lambda_1$','interpreter','latex','fontsize',15);
%% plot
for i_mode = 1:5

lon_plot = [lon(n_lon/2+1:n_lon) - 360;lon(1:n_lon/2)];
true_plot = [mode(n_lon/2+1:n_lon,:,i_mode); mode(1:n_lon/2,:,i_mode)];

load coastlines;
levels = linspace(-0.01,0.01,21);

figure;
% contourf(lon_plot,lat,true_plot',20,'linestyle','none');
% hold on;
% plot(coastlon,coastlat,'k');
ax = worldmap('world');
contourfm(lat,lon_plot,true_plot',20,'linestyle','none');
hold on;
geoshow(coastlat,coastlon,'color','k');
colorbar;
colormap(redblue);
mlabel('off');
caxis([-1e-2 1e-2]);
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
    'Position',[1 1 6 5],'PaperPosition',[0 0 9 3.5]);
outname = strcat('figs/glbal_EOF_',num2str(i_mode,'%3d'));
print(outname,'-djpeg','-r300');
close;
end