%% Figure with two panels: Lognormal PDFs and Burstiness vs BTI


% Fixed lognormal location parameter (ln X ~ N(mu, sigma^2))
mu = 0;   % adjust as you like
N=26;
% ~20 sigmas (spread of the lognormal)
sigma_vals = logspace(-1,1, N);   % >0
nSig       = numel(sigma_vals);

%% Colormap and mapping based on sigma (linear scale)
cm        = parula(N);    % base colormap
nCM       = size(cm, 1);
sigma_min = log(min(sigma_vals));
sigma_max = log(max(sigma_vals));

% Map sigma -> [0,1]
sigma_to_t = @(s) (log10(s) - sigma_min) ./ (sigma_max - sigma_min);
sigma_to_idx = @(s) max(1, min(nCM, 1 + round(sigma_to_t(s) * (nCM - 1))));

%% Precompute burstiness & BTI for each sigma
B_vals    = NaN(size(sigma_vals));
BTI_vals  = NaN(size(sigma_vals));

for k = 1:nSig
    sigma = sigma_vals(k);

    % Lognormal distribution object
    dist = makedist('Lognormal', 'mu', mu, 'sigma', sigma);

    % Coefficient of variation for lognormal: CV = std/mean
    m  = mean(dist);
    sd = std(dist);
    cv = sd ./ m;

    if ~isnan(cv) && cv > 0
        B_vals(k) = (cv - 1) ./ (cv + 1);   % burstiness in [-1,1]
    end

    % Dummy BTI function (replace with your implementation)
    BTI_vals(k) = get_BTI_from_icdf(@(x)logninv(x,mu,sigma));
end

% Mask valid points for panel B
valid        = isfinite(B_vals) & isfinite(BTI_vals);
B_plot       = B_vals(valid);
BTI_plot     = BTI_vals(valid);
sigma_plot   = sigma_vals(valid);

figure(33)
clf

%% Set up figure and layout
t = tiledlayout("horizontal",'TileSpacing', 'compact', 'Padding', 'compact');
%% Panel A: Lognormal PDFs colored by sigma

% x-grid for PDF (positive support, wide log range)
x = logspace(0, 4, 400);   % 
t1 = nexttile(t);

hold on;
for k = 1:nSig
    sigma = sigma_vals(k);
    dist  = makedist('Lognormal', 'mu', mu, 'sigma', sigma);

    y = 1-cdf(dist, x);

    % Color based on sigma
    thisColor = cm(k, :);

    plot(t1, x, y, 'LineWidth', 1.2, 'Color', thisColor);
end
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('x');
ylabel('P[X>x]');
title('(A)');
xlim([1 1E4])
ylim([1E-4 0.5])
box on;
grid(t1,'on');
%% Panel B: Burstiness vs BTI, colored by sigma
t2 = nexttile(t);

scatter(t2, B_plot, BTI_plot, 60, log10(sigma_plot), 'filled');
xlabel('BP');
ylabel('BTI');
title('(B)');
grid(t2,'on');
box on;

% Apply colormap and color limits for t in [0,1]
colormap(cm);
set(gca,'ColorScale','linear')
%clim([0 1]);

% Colorbar to the right of panel B
cb = colorbar;
cb.Label.String = 'Log_{10}(\sigma)';
cb.Location     = 'eastoutside';

fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
%t.TileSpacing = 'compact';
%t.Padding     = 'compact';
