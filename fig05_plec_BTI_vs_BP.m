%% Figure with two panels: PDFs and Burstiness vs BTI
figure(717)
clf
xmin  = 1;
alpha = 2.1;

% ~20 lambdas in logspace from 10^-4 to 10^0
lambda_vals = logspace(-4, 1, 20);
nLam        = numel(lambda_vals);

%% Common colormap and mapping based on log10(lambda)
cm      = parula(256);  % base colormap
nCM     = size(cm, 1);
lam_min = min(lambda_vals);
lam_max = max(lambda_vals);

% Map lambda -> [0,1] via log-space
lambda_to_t = @(lam) (log10(lam) - log10(lam_min)) ./ ...
                    (log10(lam_max) - log10(lam_min));
lambda_to_idx = @(lam) max(1, min(nCM, 1 + round(lambda_to_t(lam) * (nCM-1))));

%% Precompute burstiness & BTI for each lambda
B_vals   = NaN(size(lambda_vals));
BTI_vals = NaN(size(lambda_vals));

for k = 1:nLam
    lambda = lambda_vals(k);
    dist   = plcutoff_dist(alpha, lambda, xmin);

    cv = dist.cv;
    if ~isnan(cv) && cv > 0
        B_vals(k) = (cv - 1) ./ (cv + 1);      % burstiness in [-1,1]
    end

    % Dummy BTI function (replace with your implementation)
    BTI_vals(k) = get_BTI_PLEC(alpha, lambda, xmin);
end

% Mask valid points for panel B
valid = isfinite(B_vals) & isfinite(BTI_vals);
B_plot   = B_vals(valid);
BTI_plot = BTI_vals(valid);
lambda_plot = lambda_vals(valid);

%% Set up figure and layout
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Panel A: PDFs colored by lambda
t1=nexttile(1);

% x-grid for PDF (xmin=1, go out a few orders of magnitude)
x = logspace(0, 3, 400);  % from 10^0 to 10^3

hold on;
for k = 1:nLam
    lambda = lambda_vals(k);
    dist   = plcutoff_dist(alpha, lambda, xmin);

    % Assuming your dist struct has a pdf method, e.g. dist.pdf(x)
    y = 1-dist.cdf(x);

    % Color based on lambda (log-scaled)
    idx = lambda_to_idx(lambda);
    thisColor = cm(idx, :);

    plot(x, y, 'LineWidth', 2, 'Color', thisColor);
end
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('x');
ylabel('P[X>x]');
title('(A)');
box on;
grid(t1,'on')


%% Panel B: Burstiness vs BTI, colored by lambda
t2=nexttile(2);

% Use t(lambda) as color coordinate so colorbar matches mapping
t_vals = lambda_to_t(lambda_plot);

scatter(B_plot, BTI_plot, 100, t_vals, 'filled');
xlabel('BP');
ylabel('BTI');
title('(B)');
box on;
grid(t2,'on')

% Apply the colormap and color limits for t in [0,1]
colormap(t2, cm);  % applies to the whole figure / layout
%caxis([0 1]);

% Colorbar to the right of panel B
cb = colorbar;
cb.Label.String = '\lambda^*';
cb.Location = 'eastoutside';

% Choose log-spaced lambda ticks on the colorbar
lambda_ticks = [1e-4 1e-3 1e-2 1e-1 1];
lambda_ticks = lambda_ticks(lambda_ticks >= lam_min & lambda_ticks <= lam_max);
t_ticks = lambda_to_t(lambda_ticks);
cb.Ticks = t_ticks;
cb.TickLabels = arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), ...
                         lambda_ticks, 'UniformOutput', false);

% Link axes visually
t.TileSpacing = 'compact';
t.Padding     = 'compact';
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
