%% given time series of training data, construct a reduced-order model:
% long-term trend (function of global mean temperature) + sigma(t)*C * colored white noise
% long-term trends are piecewise linear
% sigma(t)^2 (variance) are also piecewise linear
% We perform a linear fitting for the mean v.s. global temperature for each season separately
% C is the cross-correlation matrix between modes
% Within each season, C is a constant matrix
% noise is generated for each season separately

% Parameters:
% data_est: data used for estimaing model parameters
% t_est: estimation / training window
% t_pred: prediction window
function data_pred = rom_glbT_lin_mean_var_pred(data_est,T_glb,T_glb_valid,t_est,t_pred)

[n_est,n_mode] = size(data_est);
n_pred = length(t_pred);
if(mod(n_est,365) > 0 || mod(n_pred,365) > 0)
    error('n_est and n_pred must be multiples of 365!!!');
    return;
end
n_year = n_est / 365;
n_year_pred = n_pred / 365;

% make sure all the vectors are column vectors
% data_est = reshape(data_est,[n_est,1]);
t_est = reshape(t_est,[n_est,1]);
t_pred = reshape(t_pred,[n_pred,1]);

%% extract seasonal mean 
% mu and T_glb are seasonally averaged and stored consecutively
% These two arrays are not used for now
coeff_season_mean = zeros(n_est,n_mode);  % not used
T_glb_season = zeros(size(T_glb));  % not used

% mu and T_glb are seasonally averaged and stored SEPARATELY
mu_est_season = cell(4,1);
mu_est_season{1} = zeros(90*n_year,n_mode); % winter
mu_est_season{2} = zeros(92*n_year,n_mode); % spring
mu_est_season{3} = zeros(92*n_year,n_mode); % summer
mu_est_season{4} = zeros(91*n_year,n_mode); % autumn
var_est_season = cell(4,1);
var_est_season{1} = zeros(90*n_year,n_mode); % winter
var_est_season{2} = zeros(92*n_year,n_mode); % spring
var_est_season{3} = zeros(92*n_year,n_mode); % summer
var_est_season{4} = zeros(91*n_year,n_mode); % autumn

T_glb_est_season = cell(4,1);
T_glb_est_season{1} = zeros(90*n_year,1); % winter
T_glb_est_season{2} = zeros(92*n_year,1); % spring
T_glb_est_season{3} = zeros(92*n_year,1); % summer
T_glb_est_season{4} = zeros(91*n_year,1); % autumn

T_glb_valid_season = cell(4,1);
T_glb_valid_season{1} = zeros(90*n_year,1); % winter
T_glb_valid_season{2} = zeros(92*n_year,1); % spring
T_glb_valid_season{3} = zeros(92*n_year,1); % summer
T_glb_valid_season{4} = zeros(91*n_year,1); % autumn

for n=1:n_year
    % winter
    if(n>1)
        i_start = (n-1)*365 - 30;
        i_end   = (n-1)*365 + 31 + 28;
        j_start = (n-1)*90 - 30;
        j_end   = (n-1)*90 + 31 + 28;
    else
        i_start = (n-1)*365 + 1;
        i_end   = (n-1)*365 + 31 + 28;
        j_start = (n-1)*90 + 1;
        j_end   = (n-1)*90 + 31 + 28;
    end
    mu_tmp = mean(data_est(i_start:i_end,:),1);
    coeff_season_mean(i_start:i_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_season(i_start:i_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    
    mu_est_season{1}(j_start:j_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_est_season{1}(j_start:j_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    T_glb_valid_season{1}(j_start:j_end) = repmat(mean(T_glb_valid(i_start:i_end)),[i_end-i_start+1,1]);

    var_tmp = var(data_est(i_start:i_end,:) - mu_tmp,[],1);
    var_est_season{1}(j_start:j_end,:) = repmat(var_tmp,[i_end-i_start+1,1]);
    
    % spring
    i_start = (n-1)*365 + 60;
    i_end   = (n-1)*365 + 151;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 31 + 30 + 31;
    mu_tmp = mean(data_est(i_start:i_end,:),1);
    coeff_season_mean(i_start:i_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_season(i_start:i_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    
    mu_est_season{2}(j_start:j_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_est_season{2}(j_start:j_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    T_glb_valid_season{2}(j_start:j_end) = repmat(mean(T_glb_valid(i_start:i_end)),[i_end-i_start+1,1]);

    var_tmp = var(data_est(i_start:i_end,:) - mu_tmp,[],1);
    var_est_season{2}(j_start:j_end,:) = repmat(var_tmp,[i_end-i_start+1,1]);
    
    % summer
    i_start = (n-1)*365 + 152;
    i_end   = (n-1)*365 + 243;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 30 + 31 + 31;
    mu_tmp = mean(data_est(i_start:i_end,:),1);
    coeff_season_mean(i_start:i_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_season(i_start:i_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    
    mu_est_season{3}(j_start:j_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_est_season{3}(j_start:j_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    T_glb_valid_season{3}(j_start:j_end) = repmat(mean(T_glb_valid(i_start:i_end)),[i_end-i_start+1,1]);

    var_tmp = var(data_est(i_start:i_end,:) - mu_tmp,[],1);
    var_est_season{3}(j_start:j_end,:) = repmat(var_tmp,[i_end-i_start+1,1]);
    
    % autumn
    i_start = (n-1)*365 + 244;
    i_end   = (n-1)*365 + 334;
    j_start = (n-1)*91 + 1;
    j_end   = (n-1)*91 + 30 + 31 + 30;
    mu_tmp = mean(data_est(i_start:i_end,:),1);
    coeff_season_mean(i_start:i_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_season(i_start:i_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    
    mu_est_season{4}(j_start:j_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
    T_glb_est_season{4}(j_start:j_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
    T_glb_valid_season{4}(j_start:j_end) = repmat(mean(T_glb_valid(i_start:i_end)),[i_end-i_start+1,1]);

    var_tmp = var(data_est(i_start:i_end,:) - mu_tmp,[],1);
    var_est_season{4}(j_start:j_end,:) = repmat(var_tmp,[i_end-i_start+1,1]);
    
    % winter in last year
    if(n==n_year)
        i_start = (n-1)*365 + 335;
        i_end   = (n-1)*365 + 365;
        j_start = (n-1)*90 + 31 + 28 + 1;
        j_end   = (n-1)*90 + 31 + 28 + 31;
        mu_tmp = mean(data_est(i_start:i_end,:),1);
        coeff_season_mean(i_start:i_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
        T_glb_season(i_start:i_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
        
        mu_est_season{1}(j_start:j_end,:) = repmat(mu_tmp,[i_end-i_start+1,1]);
        T_glb_est_season{1}(j_start:j_end) = repmat(mean(T_glb(i_start:i_end)),[i_end-i_start+1,1]);
        T_glb_valid_season{1}(j_start:j_end) = repmat(mean(T_glb_valid(i_start:i_end)),[i_end-i_start+1,1]);

        var_tmp = var(data_est(i_start:i_end,:) - mu_tmp,[],1);
        var_est_season{1}(j_start:j_end,:) = repmat(var_tmp,[i_end-i_start+1,1]);
    end
end
%% linear regression to get mu(t) = p_0 + p_1 *t
% p0 = zeros(n_mode,1);
% p1 = zeros(n_mode,1);
% data_est_fluct = zeros(size(data_est));
% % for j=1:n_mode
% %     mdl = polyfit(T_glb,data_est(:,j),1);
% %     p1(j) = mdl(1);
% %     p0(j) = mdl(2);
% %     data_est_fluct(:,j) = data_est(:,j) - (p0(j) + p1(j)*T_glb);
% % end

% for j=1:n_mode
%     mdl = polyfit(T_glb_season,coeff_season_mean(:,j),1);
%     p1(j) = mdl(1);
%     p0(j) = mdl(2);
%     data_est_fluct(:,j) = data_est(:,j) - (p0(j) + p1(j)*T_glb_season);
% end

p0_mu = zeros(n_mode,4);
p1_mu = zeros(n_mode,4);
p0_var = zeros(n_mode,4);
p1_var = zeros(n_mode,4);
% data_est_fluct = zeros(size(data_est));

mu_pred_season = cell(4,1);
mu_pred_season{1} = zeros(90*n_year,n_mode); % winter
mu_pred_season{2} = zeros(92*n_year,n_mode); % spring
mu_pred_season{3} = zeros(92*n_year,n_mode); % summer
mu_pred_season{4} = zeros(91*n_year,n_mode); % autumn

mu_valid_season = cell(4,1);
mu_valid_season{1} = zeros(90*n_year,n_mode); % winter
mu_valid_season{2} = zeros(92*n_year,n_mode); % spring
mu_valid_season{3} = zeros(92*n_year,n_mode); % summer
mu_valid_season{4} = zeros(91*n_year,n_mode); % autumn

var_pred_season = cell(4,1);
var_pred_season{1} = zeros(90*n_year,n_mode); % winter
var_pred_season{2} = zeros(92*n_year,n_mode); % spring
var_pred_season{3} = zeros(92*n_year,n_mode); % summer
var_pred_season{4} = zeros(91*n_year,n_mode); % autumn

var_valid_season = cell(4,1);
var_valid_season{1} = zeros(90*n_year,n_mode); % winter
var_valid_season{2} = zeros(92*n_year,n_mode); % spring
var_valid_season{3} = zeros(92*n_year,n_mode); % summer
var_valid_season{4} = zeros(91*n_year,n_mode); % autumn

for k=1:4
    for j=1:n_mode
        mdl = polyfit(T_glb_est_season{k},mu_est_season{k}(:,j),1);
        p1_mu(j,k) = mdl(1);
        p0_mu(j,k) = mdl(2);
        mu_pred_season{k}(:,j) = (p0_mu(j,k) + p1_mu(j,k)*T_glb_est_season{k});
        mu_valid_season{k}(:,j) = (p0_mu(j,k) + p1_mu(j,k)*T_glb_valid_season{k});

        mdl = polyfit(T_glb_est_season{k},var_est_season{k}(:,j),1);
        p1_var(j,k) = mdl(1);
        p0_var(j,k) = mdl(2);
        var_pred_season{k}(:,j) = (p0_var(j,k) + p1_var(j,k)*T_glb_est_season{k});
        var_valid_season{k}(:,j) = (p0_var(j,k) + p1_var(j,k)*T_glb_valid_season{k});
    end
end

mu_pred = zeros(n_est,n_mode);
var_pred = zeros(n_est,n_mode);
mu_valid = zeros(n_est,n_mode);
var_valid = zeros(n_est,n_mode);

for n=1:n_year
    % winter
    i_start = (n-1)*365 + 1;
    i_end   = (n-1)*365 + 31 + 28;
    j_start = (n-1)*90 + 1;
    j_end   = (n-1)*90 + 31 + 28;
    mu_pred(i_start:i_end,:) = mu_pred_season{1}(j_start:j_end,:);
    var_pred(i_start:i_end,:) = var_pred_season{1}(j_start:j_end,:);
    mu_valid(i_start:i_end,:) = mu_valid_season{1}(j_start:j_end,:);
    var_valid(i_start:i_end,:) = var_valid_season{1}(j_start:j_end,:);
    % spring
    i_start = (n-1)*365 + 60;
    i_end   = (n-1)*365 + 151;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 31 + 30 + 31;
    mu_pred(i_start:i_end,:) = mu_pred_season{2}(j_start:j_end,:);
    var_pred(i_start:i_end,:) = var_pred_season{2}(j_start:j_end,:);
    mu_valid(i_start:i_end,:) = mu_valid_season{2}(j_start:j_end,:);
    var_valid(i_start:i_end,:) = var_valid_season{2}(j_start:j_end,:);
    % summer
    i_start = (n-1)*365 + 152;
    i_end   = (n-1)*365 + 243;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 30 + 31 + 31;
    mu_pred(i_start:i_end,:) = mu_pred_season{3}(j_start:j_end,:);
    var_pred(i_start:i_end,:) = var_pred_season{3}(j_start:j_end,:);
    mu_valid(i_start:i_end,:) = mu_valid_season{3}(j_start:j_end,:);
    var_valid(i_start:i_end,:) = var_valid_season{3}(j_start:j_end,:);
    % autumn
    i_start = (n-1)*365 + 244;
    i_end   = (n-1)*365 + 334;
    j_start = (n-1)*91 + 1;
    j_end   = (n-1)*91 + 30 + 31 + 30;
    mu_pred(i_start:i_end,:) = mu_pred_season{4}(j_start:j_end,:);
    var_pred(i_start:i_end,:) = var_pred_season{4}(j_start:j_end,:);
    mu_valid(i_start:i_end,:) = mu_valid_season{4}(j_start:j_end,:);
    var_valid(i_start:i_end,:) = var_valid_season{4}(j_start:j_end,:);
    % winter
    i_start = (n-1)*365 + 335;
    i_end   = (n-1)*365 + 365;
    j_start = (n-1)*90 + 31 + 28 + 1;
    j_end   = (n-1)*90 + 31 + 28 + 31;
    mu_pred(i_start:i_end,:) = mu_pred_season{1}(j_start:j_end,:);
    var_pred(i_start:i_end,:) = var_pred_season{1}(j_start:j_end,:);
    mu_valid(i_start:i_end,:) = mu_valid_season{1}(j_start:j_end,:);
    var_valid(i_start:i_end,:) = var_valid_season{1}(j_start:j_end,:);
end

% check if there is any negative variance
idx = find(var_pred(:) < 0);
[idx_1,idx_2] = ind2sub(size(var_pred),idx);
if(~isempty(idx))
    fprintf('negative variance found at %d, %d\n',idx_1,idx_2);
    return;
end
data_est_fluct = (data_est - mu_pred) ./ sqrt(var_pred);
% data_est_fluct = (data_est - coeff_season_mean) ./ sqrt(var_pred); % use true mean

%% separate data for each season
data_est_season = cell(4,1);
data_est_season{1} = zeros(90*n_year,n_mode); % winter
data_est_season{2} = zeros(92*n_year,n_mode); % spring
data_est_season{3} = zeros(92*n_year,n_mode); % summer
data_est_season{4} = zeros(91*n_year,n_mode); % autumn
for n=1:n_year
    i_start = (n-1)*365 + 1;
    i_end   = (n-1)*365 + 31 + 28;
    j_start = (n-1)*90 + 1;
    j_end   = (n-1)*90 + 31 + 28;
    data_est_season{1}(j_start:j_end,:) = data_est_fluct(i_start:i_end,:);
    i_start = (n-1)*365 + 60;
    i_end   = (n-1)*365 + 151;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 31 + 30 + 31;
    data_est_season{2}(j_start:j_end,:) = data_est_fluct(i_start:i_end,:);
    i_start = (n-1)*365 + 152;
    i_end   = (n-1)*365 + 243;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 30 + 31 + 31;
    data_est_season{3}(j_start:j_end,:) = data_est_fluct(i_start:i_end,:);
    i_start = (n-1)*365 + 244;
    i_end   = (n-1)*365 + 334;
    j_start = (n-1)*91 + 1;
    j_end   = (n-1)*91 + 30 + 31 + 30;
    data_est_season{4}(j_start:j_end,:) = data_est_fluct(i_start:i_end,:);
    i_start = (n-1)*365 + 335;
    i_end   = (n-1)*365 + 365;
    j_start = (n-1)*90 + 31 + 28 + 1;
    j_end   = (n-1)*90 + 31 + 28 + 31;
    data_est_season{1}(j_start:j_end,:) = data_est_fluct(i_start:i_end,:);
end
% compute correlation and Cholesky decomposition (Eq.(9) in the NeurIPS paper)
corr_season = zeros(n_mode,n_mode,4);
L_season = zeros(n_mode,n_mode,4);  % L_season is C_s in Eq.(9)
L_inv_season = zeros(n_mode,n_mode,4); % C_s^{-1}
for i=1:4
    [l1,~] = size(data_est_season{i});
    corr_season(:,:,i) = data_est_season{i}' * data_est_season{i} / l1;
    L_season(:,:,i) = chol(corr_season(:,:,i))';
    L_inv_season(:,:,i) = inv(L_season(:,:,i));
end
% extract the "noise" part by transforming the data with L_inv
data_est_noise = cell(4,1);
data_est_noise{1} = zeros(90*n_year,n_mode); % winter
data_est_noise{2} = zeros(92*n_year,n_mode); % spring
data_est_noise{3} = zeros(92*n_year,n_mode); % summer
data_est_noise{4} = zeros(91*n_year,n_mode); % autumn
% Eq.(10) in the NeurIPS paper
for i=1:4
    tmp = L_inv_season(:,:,i) * data_est_season{i}';
    data_est_noise{i} = tmp';  % w_s(t) in Eq. (10)
end
clear data_est_season;
clear data_est_fluct;

% compute the spectrum of "colored noise" and simulate the noise
ratio = floor(n_pred/n_est);
residual = n_pred - n_est * ratio;
% data_pred_noise = zeros(n_pred,n_mode);

data_pred_season = cell(4,1);
data_pred_season{1} = zeros(90*n_year,n_mode); % winter
data_pred_season{2} = zeros(92*n_year,n_mode); % spring
data_pred_season{3} = zeros(92*n_year,n_mode); % summer
data_pred_season{4} = zeros(91*n_year,n_mode); % autumn

n_length = [90,92,92,91]*n_year;

data_pred_fluct = zeros(n_pred,n_mode);

%% generate colored noise
season_days = [90,92,92,91];
for k=1:4   
    for j_mode = 1:n_mode
        % compute spectrum
        if(k==1)
            spectrum = zeros(season_days(k),1);
            for n=1:n_year-1
                i_start = (n-1)*season_days(k) + 31 + 28 + 1;
                i_end = n*season_days(k) + 31 + 28;
                spectrum = spectrum + abs(fft(data_est_noise{k}(i_start:i_end,j_mode))).^2;
            end
            spectrum = spectrum / n_year;
        else
            spectrum = zeros(season_days(k),1);
            for n=1:n_year
                i_start = (n-1)*season_days(k) + 1;
                i_end = n*season_days(k);
                spectrum = spectrum + abs(fft(data_est_noise{k}(i_start:i_end,j_mode))).^2;
            end
            spectrum = spectrum / n_year;
        end
        % given autocorrelation, generate samples
        autocorr = ifft(spectrum) / season_days(k);
        if(k==1)
            for n=1:n_year-1
                i_start = (n-1)*season_days(k) + 31 + 28 + 1;
                i_end = n*season_days(k) + 31 + 28;
                Y_t = autocorr_to_colored_noise(autocorr,season_days(k));
                data_pred_season{k}(i_start:i_end,j_mode) = Y_t;
            end
            Y_t = autocorr_to_colored_noise(autocorr,season_days(k));
            data_pred_season{k}(1:59,j_mode) = Y_t(1:59);
            data_pred_season{k}(end-30:end,j_mode) = Y_t(60:90);
        else
            for n=1:n_year_pred
                i_start = (n-1)*season_days(k) + 1;
                i_end = n*season_days(k);
                Y_t = autocorr_to_colored_noise(autocorr,season_days(k));
                data_pred_season{k}(i_start:i_end,j_mode) = Y_t;
            end
        end

    end
    % Multiply the colored noise by L_season
    tmp = L_season(:,:,k) * data_pred_season{k}';
    data_pred_season{k} = tmp';
end

offset = 0;

for n=1:n_year
    % winter
    i_start = (n-1)*365 + 1 + offset;
    i_end   = (n-1)*365 + 31 + 28 + offset;
    j_start = (n-1)*90 + 1;
    j_end   = (n-1)*90 + 31 + 28;
    data_pred_fluct(i_start:i_end,:) = data_pred_season{1}(j_start:j_end,:);
    % spring
    i_start = (n-1)*365 + 60 + offset;
    i_end   = (n-1)*365 + 151 + offset;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 31 + 30 + 31;
    data_pred_fluct(i_start:i_end,:) = data_pred_season{2}(j_start:j_end,:);
    % summer
    i_start = (n-1)*365 + 152 + offset;
    i_end   = (n-1)*365 + 243 + offset;
    j_start = (n-1)*92 + 1;
    j_end   = (n-1)*92 + 30 + 31 + 31;
    data_pred_fluct(i_start:i_end,:) = data_pred_season{3}(j_start:j_end,:);
    % autumn
    i_start = (n-1)*365 + 244 + offset;
    i_end   = (n-1)*365 + 334 + offset;
    j_start = (n-1)*91 + 1;
    j_end   = (n-1)*91 + 30 + 31 + 30;
    data_pred_fluct(i_start:i_end,:) = data_pred_season{4}(j_start:j_end,:);
    % winter
    i_start = (n-1)*365 + 335 + offset;
    i_end   = (n-1)*365 + 365 + offset;
    j_start = (n-1)*90 + 31 + 28 + 1;
    j_end   = (n-1)*90 + 31 + 28 + 31;
    data_pred_fluct(i_start:i_end,:) = data_pred_season{1}(j_start:j_end,:);
end



clear data_est_noise;

% add back linear part
% data_pred = zeros(n_pred,n_mode);
% for j_mode = 1:n_mode
%     data_pred(:,j_mode) = data_pred_fluct(:,j_mode) + ...
%         p0(j_mode) + p1(j_mode) * T_glb_season(t_pred); % Eq.(11)
% end
data_pred = data_pred_fluct.*sqrt(var_valid) + mu_valid;
% data_pred = data_pred_fluct.*sqrt(var_pred) + coeff_season_mean;




