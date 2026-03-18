%% Heat map of burstiness over (alpha, lambda) with log-spaced lambda
clear; clc;

% Parameter grid
alpha_vals  = linspace(0, 4, 200)+1E-4; % unstable solution for integer values, so add +1E4
lambda_vals = logspace(-4, 0, 200); % log-spaced lambda from 1e-3 to 1e1
xmin = 1;

[Agrid, Lgrid] = meshgrid(alpha_vals, lambda_vals);
Bmap = NaN(size(Agrid));

%% Compute burstiness for each pair (alpha, lambda)
for i = 1:length(lambda_vals)
    for j = 1:length(alpha_vals)
        alpha  = Agrid(i,j);
        lambda = Lgrid(i,j);

        dist = plcutoff_dist(alpha, lambda, xmin);
        cv   = dist.cv;

        if isnan(cv) || cv <= 0
            Bmap(i,j) = NaN;
        else
            % burstiness measure
            Bmap(i,j) = (cv - 1) ./ (cv + 1);   % range [-1, 1]
        end
    end
end

%% Diverging custom colormap (SWAPPED: red for negative, blue for positive)
n = 256;

% negative values colormap --> RED
neg_cm = ([linspace(0.7,1,n/2)' ...   % R ↑
                 linspace(0.2,0.9,n/2)' ... % G ↑
                 linspace(0,0.8,n/2)']);    % B ↑

% positive values colormap --> BLUE
pos_cm = flipud([linspace(0,0.8,n/2)' ...          % R ↑
          linspace(0.2,0.9,n/2)' ...        % G ↑
          linspace(0.7,1,n/2)']);            % B ↑

cmap = [neg_cm; pos_cm];

%% Plot heat map
figure;
imagesc(alpha_vals, lambda_vals, Bmap);
set(gca,'YDir','normal');

set(gca,'YScale','log');      % log λ axis
colormap(cmap);
caxis([-1 1]);

c = colorbar;
c.Label.String = 'BP';

xlabel('\alpha');
ylabel('\lambda');
title('');

% Zero contour
hold on;
contour(alpha_vals, lambda_vals, Bmap, [0 0], 'k-', 'LineWidth', 1.2);
hold off;
