%% modified from coeff_global_spectrum.m
clear;
close all;

n_days_season = 92;
n_rlz = 10;

load('pca_historical/modes_Tmx_365skip5_global_historical.mat','lambda');

load('true_rom_w.mat','true_w');
true_w_summer = true_w{3};
clear true_w;
load('rom_w_summer.mat');

% WARNING!!!!!!!!!!!!!!!!
% % use mode 1, 2, 500
j_mode = 500;
true_w = zeros(251*n_days_season,3);
true_w(:,1:2) = true_w_summer(:,1:2);
true_w(:,3) = true_w_summer(:,j_mode);
rom_w_mrlz = zeros(251*n_days_season,3,n_rlz);
rom_w_mrlz(:,1:2,:) = rom_w_summer_mrlz(:,1:2,:);
rom_w_mrlz(:,3,:) = rom_w_summer_mrlz(:,j_mode,:);
clear true_w_summer rom_w_summer_mrlz;
%% historical spectrum
true_w2_f_hist = zeros(n_days_season,3);
rom_w2_f_hist = zeros(n_days_season,3);

kstart = 1;  %170;
kend = 165;  %250;
rom_w2_f_mrlz = zeros(n_days_season,3,n_rlz);
for k=kstart:kend
    jstart = (k-1)*n_days_season+1;
    jend   = k*n_days_season;
    tmp = fft(true_w(jstart:jend,:),[],1)/n_days_season;
    true_w2_f_hist = true_w2_f_hist + abs(tmp).^2;

    for l=1:n_rlz
        tmp = fft(rom_w_mrlz(jstart:jend,:,l),[],1)/n_days_season;
        rom_w2_f_mrlz(:,:,l) = rom_w2_f_mrlz(:,:,l) + abs(tmp).^2;
        rom_w2_f_hist = rom_w2_f_hist + abs(tmp).^2;
    end
end
true_w2_f_hist = true_w2_f_hist / (kend - kstart + 1);
rom_w2_f_hist = rom_w2_f_hist / (kend - kstart + 1) / n_rlz;

%% ssp spectrum
true_w2_f_ssp = zeros(n_days_season,3);
rom_w2_f_ssp = zeros(n_days_season,3);

kstart = 166;  %170;
kend = 251;  %250;
rom_w2_f_mrlz = zeros(n_days_season,3,n_rlz);
for k=kstart:kend
    jstart = (k-1)*n_days_season+1;
    jend   = k*n_days_season;
    tmp = fft(true_w(jstart:jend,:),[],1)/n_days_season;
    true_w2_f_ssp = true_w2_f_ssp + abs(tmp).^2;

    for l=1:n_rlz
        tmp = fft(rom_w_mrlz(jstart:jend,:,l),[],1)/n_days_season;
        rom_w2_f_mrlz(:,:,l) = rom_w2_f_mrlz(:,:,l) + abs(tmp).^2;
        rom_w2_f_ssp = rom_w2_f_ssp + abs(tmp).^2;
    end
end
true_w2_f_ssp = true_w2_f_ssp / (kend - kstart + 1);
rom_w2_f_ssp = rom_w2_f_ssp / (kend - kstart + 1) / n_rlz;

%% only keep the positive wavenumber componets
L = n_days_season;
epsilon = 1e-10;
true_w2_f_hist = true_w2_f_hist(1:L/2+1,:);
true_w2_f_hist(2:end-1,:) = 2*true_w2_f_hist(2:end-1,:);
rom_w2_f_hist = rom_w2_f_hist(1:L/2+1,:);
rom_w2_f_hist(2:end-1,:) = 2*rom_w2_f_hist(2:end-1,:);
true_w2_f_ssp = true_w2_f_ssp(1:L/2+1,:);
true_w2_f_ssp(2:end-1,:) = 2*true_w2_f_ssp(2:end-1,:);
rom_w2_f_ssp = rom_w2_f_ssp(1:L/2+1,:);
rom_w2_f_ssp(2:end-1,:) = 2*rom_w2_f_ssp(2:end-1,:);

%% plot
frequency = (1:(L/2-1))/L;
period = 1 ./ frequency;

figure;

% first panel
i_mode = 1;
subplot(1,3,1);
plot(period,true_w2_f_hist(2:end-1,i_mode),'o','Color',[0.5 0.5 0.5],'linewidth',1.5);
hold on;
plot(period,true_w2_f_ssp(2:end-1,i_mode),'Color',[1,0,0],'linewidth',2);
plot(period,rom_w2_f_ssp(2:end-1,i_mode),...
    'Color',[0,0,1,0.8],'linewidth',1.8);
set(gca,'fontsize',13);
set(gca,'TickLabelInterpreter', 'latex','YScale','Log','XScale','Log');
xlabel('days','fontsize',15,'fontname','Times');
ylabel('$|\tilde{\eta}_i|^2$','interpreter','latex','fontsize',15);
axis([1 1e2 1e-4 1e-0]);
text(-0.25,0.95,'($a$)','interpreter','latex','Units','normalized','FontSize',14);
h = legend('Historical data','SSP585 data','Emulation','location','northwest');
set(h,'box','off','fontsize',9);
% set(h,'interpreter','latex');

% second panel
i_mode = 2;
subplot(1,3,2);
plot(period,true_w2_f_hist(2:end-1,i_mode),'o','Color',[0.5 0.5 0.5],'linewidth',1.5);
hold on;
plot(period,true_w2_f_ssp(2:end-1,i_mode),'Color',[1,0,0],'linewidth',2);
plot(period,rom_w2_f_ssp(2:end-1,i_mode),...
    'Color',[0,0,1,0.8],'linewidth',1.8);
set(gca,'fontsize',13);
set(gca,'TickLabelInterpreter', 'latex','YScale','Log','XScale','Log');
xlabel('days','fontsize',15,'fontname','Times');
ylabel('$|\tilde{\eta}_i|^2$','interpreter','latex','fontsize',15);
axis([1 1e2 1e-4 1e-0]);
text(-0.25,0.95,'($b$)','interpreter','latex','Units','normalized','FontSize',14);

% third panel
i_mode = 3;
subplot(1,3,3);
plot(period,true_w2_f_hist(2:end-1,i_mode),'o','Color',[0.5 0.5 0.5],'linewidth',1.5);
hold on;
plot(period,true_w2_f_ssp(2:end-1,i_mode),'Color',[1,0,0],'linewidth',2);
plot(period,rom_w2_f_ssp(2:end-1,i_mode),...
    'Color',[0,0,1,0.8],'linewidth',1.8);
set(gca,'fontsize',13);
set(gca,'TickLabelInterpreter', 'latex','YScale','Log','XScale','Log');
xlabel('days','fontsize',15,'fontname','Times');
ylabel('$|\tilde{\eta}_i|^2$','interpreter','latex','fontsize',15);
axis([1 1e2 1e-4 1e-0]);
text(-0.25,0.95,'($c$)','interpreter','latex','Units','normalized','FontSize',14);

set(gcf,'Units','inches',...
        'Position',[1 1 13 2.3]);
