% DEPENDENCIES: randht by Aaron Clauset - 
% https://aaronclauset.github.io/powerlaws/


% Plot MSE surfaces for BP/BTI under time-window censoring.
% Generates a 2x2 subplot figure:
%   {power law, power law w/ cutoff, lognormal, Weibull} with BP/BTI overlaid
%
% NOTE: This script assumes xmin = 1 for the power-law family by default.

clear; 
figure(723)
clc;

%% Experiment controls
Npts = 25;
N_samples = 1E3;   % number of samples per trial
M_trials = 1E3;    % number of trials
T_WINDOW = 1E3;     % time-window censoring parameter
xmin = 1;          % power-law minimum value (randht default)

%% Parameter ranges (edit as needed)
% Power law (single parameter: alpha)
alpha_pl_vals = linspace(2.0001,2+sqrt(2),Npts);

% Power law w/ cutoff (PLEC): fix alpha, vary lambda
alpha_cut = 2.5;
lambda_cut_vals = logspace(-4,-1.3,Npts);

% Lognormal (two parameters: mu, sigma)
mu_ln = 1;
sigma_ln_vals = linspace(sqrt(log(2)),5,Npts);

% Weibull (two parameters: scale lambda, shape k)
lambda_wb = 1;
k_wb_vals = linspace(0.1,2.5,Npts);

%% Allocate MSE arrays
mse_bp_pl = zeros(size(alpha_pl_vals));
mse_bti_pl = zeros(size(alpha_pl_vals));

mse_bp_cut = zeros(size(lambda_cut_vals));
mse_bti_cut = zeros(size(lambda_cut_vals));

mse_bp_ln = zeros(size(sigma_ln_vals));
mse_bti_ln = zeros(size(sigma_ln_vals));

mse_bp_wb = zeros(size(k_wb_vals));
mse_bti_wb = zeros(size(k_wb_vals));

%% Helper: time-window censoring
censor_samples = @(x) min(x, T_WINDOW * rand(size(x)));

%% Power law (Pareto Type I, xmin = 1)
for ii = 1:numel(alpha_pl_vals)
    alpha = alpha_pl_vals(ii);

    % analytical mean/variance for Clauset-style power law:
    % f(x) = (alpha-1) * xmin^(alpha-1) * x^(-alpha)
    mu = ((alpha - 1) * xmin) / (alpha - 2);
    varx = ((alpha - 1) * xmin^2) / ((alpha - 3) * (alpha - 2)^2);
    r = sqrt(varx) / mu;
    BP = (r - 1) / (r + 1);
    BP(alpha<3)=1;

    icdf_pl = @(p) xmin .* (1 - p).^(-1 ./ (alpha - 1));
    BTI = get_BTI_from_icdf(icdf_pl);

    x = randht(M_trials * N_samples, 'powerlaw', alpha, 'xmin', xmin);
    x = reshape(x, M_trials, N_samples);
    y = censor_samples(x);

    [mse_bp, mse_bti] = mse_sample_size_fresh_samples(y, BP, BTI);
    mse_bp_pl(ii) = mse_bp;
    mse_bti_pl(ii) = mse_bti;
end

%% Power law with cutoff (PLEC)
for jj = 1:numel(lambda_cut_vals)
    lambda = lambda_cut_vals(jj);

    BP = get_PLEC_BP(alpha_cut, lambda, xmin);
    BTI = get_BTI_PLEC(alpha_cut, lambda, xmin);

    x = randht(M_trials * N_samples, 'cutoff', alpha_cut, lambda, 'xmin', xmin);
    x = reshape(x, M_trials, N_samples);
    y = censor_samples(x);

    [mse_bp, mse_bti] = mse_sample_size_fresh_samples(y, BP, BTI);
    mse_bp_cut(jj) = mse_bp;
    mse_bti_cut(jj) = mse_bti;
end

%% Lognormal
for jj = 1:numel(sigma_ln_vals)
    sigma = sigma_ln_vals(jj);

    [mu_x, var_x] = lognstat(mu_ln, sigma);
    r = sqrt(var_x) / mu_x;
    BP = (r - 1) / (r + 1);

    icdf_ln = @(p) logninv(p, mu_ln, sigma);
    BTI = get_BTI_from_icdf(icdf_ln);

    x = lognrnd(mu_ln, sigma, M_trials, N_samples);
    y = censor_samples(x);

    [mse_bp, mse_bti] = mse_sample_size_fresh_samples(y, BP, BTI);
    mse_bp_ln(jj) = mse_bp;
    mse_bti_ln(jj) = mse_bti;
end

%% Weibull (stretched exponential)
for jj = 1:numel(k_wb_vals)
    k = k_wb_vals(jj);

    [mu_x, var_x] = wblstat(lambda_wb, k);
    r = sqrt(var_x) / mu_x;
    BP = (r - 1) / (r + 1);

    icdf_wb = @(p) wblinv(p, lambda_wb, k);
    BTI = get_BTI_from_icdf(icdf_wb);

    x = wblrnd(lambda_wb, k, M_trials, N_samples);
    y = censor_samples(x);

    [mse_bp, mse_bti] = mse_sample_size_fresh_samples(y, BP, BTI);
    mse_bp_wb(jj) = mse_bp;
    mse_bti_wb(jj) = mse_bti;
end

%% Plotting

%figure(723, 'Color', 'w');
layout = tiledlayout(1, 3, 'TileSpacing', 'compact');

% Power law
nexttile(layout, 1);
plot(alpha_pl_vals, sqrt(mse_bp_pl), '-o', ...
    alpha_pl_vals, sqrt(mse_bti_pl), ':+', ...
    'LineWidth', 3);
xlabel('\alpha'); ylabel('RMSE');
title('Power law');
grid on
legend('BP', 'BTI', 'Location', 'best');

% Power law w/ cutoff
nexttile(layout, 2);
semilogx(lambda_cut_vals, sqrt(mse_bp_cut), '-o', ...
    lambda_cut_vals, sqrt(mse_bti_cut), ':+', ...
    'LineWidth', 3);
xlabel('\lambda'); ylabel('RMSE');
title(sprintf('Power law w/ cutoff (\\alpha = %.2f)', alpha_cut));
grid on
legend('BP', 'BTI', 'Location', 'best');

% Lognormal
nexttile(layout, 3);
plot(sigma_ln_vals, sqrt(mse_bp_ln), '-o', ...
    sigma_ln_vals, sqrt(mse_bti_ln), ':+', ...
    'LineWidth', 3);
xlabel('\sigma'); ylabel('RMSE');
title(sprintf('Lognormal (\\mu = %.2f)', mu_ln));
grid on
legend('BP', 'BTI', 'Location', 'best');


% Weibull
% nexttile(layout, 4);
% plot(k_wb_vals, sqrt(mse_bp_wb), '-o', k_wb_vals, sqrt(mse_bti_wb), ':+', ...
%     'LineWidth', 3);
% xlabel('k'); ylabel('RMSE');
% title(sprintf('Weibull (stretched, \\lambda = %.2f)', lambda_wb));
% grid on
% legend('BP', 'BTI', 'Location', 'best');
