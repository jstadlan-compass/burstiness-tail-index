% Written with chatgpt
% Define the range for x values
phi = 0.99;
plo = 0.75;

alpha_low = 1;
alpha_hi = 4;

x = linspace(1, 10, 1000);  % Pareto distribution with x_min = 1

zero_burstiness = 2+sqrt(2);

% Define the alpha values
alpha_values = [alpha_low,2,1+sqrt(2),3,alpha_hi];

% Initialize the figure
figure(19);
clf;

% Left subplot: Complementary CDF (1 - CDF)
subplot(1, 2, 1);
hold on;

for alpha = alpha_values
    % Compute the complementary CDF
    cdf_values = 1 - (1 ./ x).^alpha;
    complementary_cdf = 1 - cdf_values;
    
    % Plot the complementary CDF
    loglog(x, complementary_cdf, 'DisplayName', "\alpha=" + num2str(alpha+1),'LineWidth',2);
end

xlabel('x');
ylabel('P[X>x]');
%title('Complementary CDF of Pareto Distribution, x_{min}=1');
legend('show');
set(gca,'XScale','log')
set(gca,'YScale','log')
grid on
hold off;

% Right subplot: Coefficient of variation and quantile ratio
subplot(1, 2, 2);
hold on;

alpha_range = linspace(alpha_low-1,alpha_hi,10000);
alpha_range_b = [2 linspace(2+eps,alpha_hi, 10000)];

% Calculate coefficient of variation and quantile ratio
cv = @(alpha) 1./sqrt(alpha.*(alpha-2));  % Coefficient of variation
q99 = @(alpha) (1 - phi).^(-1./alpha);  % 0.99 quantile
q75 = @(alpha) (1 - plo).^(-1./alpha);  % 0.75 quantile
qmed = @(alpha) (1 - .5).^(-1./alpha);  % 0.75 quantile

cv_values = cv(alpha_range);
burstiness = (cv_values-1)./(cv_values+1);

quantile_ratio = (q99(alpha_range)-qmed(alpha_range)) ./ (q75(alpha_range)-qmed(alpha_range));

quantile_exp_ref = (-log(1-phi)-log(2))/(-log(1-plo)-log(2));
r = quantile_ratio / quantile_exp_ref;
tail_ratio = (r-1)./(r+1);

% Plot the coefficient of variation
ind_burst_valid = alpha_range >= 2;
ylims = [min(burstiness(ind_burst_valid)) max(tail_ratio)+.2];
plot([2 3 alpha_range(ind_burst_valid)+1], [1 1 burstiness(ind_burst_valid)], 'r', 'DisplayName', 'BP','LineWidth',3);

% Plot the quantile ratio
plot(alpha_range+1, tail_ratio, 'b', 'DisplayName', 'BTI','LineWidth',3);
plot(2+sqrt(2),0,'rs', 'DisplayName', 'BP(2+2^{1/2})','LineWidth',4);
%plot([zero_burstiness zero_burstiness],ylims,'r--', 'DisplayName', 'BP(x) = 0','LineWidth',2);
plot([1 5],[-1 -1],'r:', 'DisplayName', 'lower limit of BP(\alpha)','LineWidth',2);
plot([1 5],[0 0],'b:', 'DisplayName', 'lower limit of BTI(\alpha)','LineWidth',2);

plot(1,1,'bo','DisplayName','','LineWidth',3);
plot(2,1,'ro','DisplayName','','LineWidth',3);

legend('BP','BTI','BP(2+2^{1/2}','lower limit of BP(\alpha)','lower limit of BTI(\alpha)','','');


%title('Burstiness Metrics (Analytical)')
grid on
xlim([alpha_low alpha_hi])

% from chatgpt:
specific_x = zero_burstiness;
specific_label = '(2+2^{1/2})';

% Get the current x-ticks and x-tick labels
%current_xticks = xticks;
%current_xticklabels = xticklabels;
current_xticks = 1:5;
current_xticklabels = string((1:5)');

% Add the specific x-value to the list of x-ticks
new_xticks = [current_xticks, specific_x];

% Add the specific label to the list of x-tick labels
new_xticklabels = [current_xticklabels; {specific_label}];

[new_xticks,ind_sorted] = sort(new_xticks);
new_xticklabels = new_xticklabels(ind_sorted);

% Set the new x-ticks and x-tick labels
%set(gca, 'TickLabelInterpreter', 'latex');
xticks(new_xticks);
xticklabels(new_xticklabels);

ylim(ylims)
xlabel('\alpha');
ylabel('Parameter Value on [-1,1]');
ylim([-0.4 1])
xlim([1 5])
hold off;
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")