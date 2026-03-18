figure(394)
clf

% Shared parameters
p_BTI = [.5 .75 .99];
rangeEnd = 1E4;
colorSet = lines(8);

doMedian = true;
if doMedian
    avgFun = @(x) median(x);
    disp("normalize by median")
else
    avgFun = @(x) mean(x);
    disp("normalize by mean")
end

%% Column 1: Sample-size plots
N_samples = 1E5;
M_samples = 1E3;

x_values_window = logspace(0, 5, N_window);
xlims_window = [1E1 1E4];
ylims_window = [-0.5 1];

x_values_samples = logspace(0, 4, N_samples);
xlims_samples = [10 N_samples];
ylims_samples = [-1 1];

subplot(3,2,1)
% Power Law
type = 'pl';
alpha_PL = 2.1;

gamma_var = alpha_PL + 1; % convert from PL to Pareto
k_gpd = 1/gamma_var;
sigma_gpd = gamma_var;
theta_gpd = 1;
gpd = makedist('GeneralizedPareto','k',k_gpd,'sigma',sigma_gpd,'theta',theta_gpd);
mean_val = mean(gpd);
std_val = std(gpd);
quants = icdf(gpd,p_BTI);
samples = random(gpd,M_samples,N_samples);
CoV = std_val / mean_val;
BP = (CoV - 1) / (CoV + 1);
BTI = get_tail_ratio_from_quantiles(quants,p_BTI);

plot_per_sample_size(samples, colorSet(1,:),colorSet(2,:),BP,BTI);
set(gca, 'XScale', 'log', 'YScale', 'linear');
xlim(xlims_samples);
ylim(ylims_samples)
legend("BP-sampled (mean)","BP-sampled (75% interval)","BP-analytical","BTI-sampled (mean)","BTI-sampled (75% interval)","BTI-analytical",'location','southeast')
title("Power Law, \alpha=2.1")
xlabel("Number of Samples")
ylabel("Burstiness Metric Value")

% Time window
samples = samples(:);
avg = avgFun(samples);
subplot(3,2,2)
disp(type + ": " + num2str(avg))
plot_windowed(samples/avg, colorSet(1,:),colorSet(2,:),BP,BTI);
set(gca, 'XScale', 'log', 'YScale', 'linear');
xlim(xlims_window);
ylim(ylims_window)
legend("BP-windowed","BTI-windowed","BP-infinite","BTI-infinite",'location','southeast')
title("Power Law, \alpha=2.1")
xlabel("Time Window T_f");
ylabel("Burstiness Metric Value");


% Power Law with Exponential Cut-Off
subplot(3,2,3)
type = 'plcut';
gamma_var = 2.5;
lambda = 1E-2;

samples = randht(M_samples * N_samples, 'cutoff', gamma_var, lambda);
samples = reshape(samples, M_samples, N_samples);
pldist = plcutoff_dist(gamma_var, lambda, 1);

mean_val = pldist.mean;
std_val = sqrt(pldist.var);
qhi = pldist.ppf(p_BTI(3));
qlo = pldist.ppf(p_BTI(2));
qmedian = pldist.ppf(p_BTI(1));
quants = [qmedian qlo qhi];
CoV = std_val / mean_val;
BP = (CoV - 1) / (CoV + 1);
BTI = get_tail_ratio_from_quantiles(quants,p_BTI);

plot_per_sample_size(samples, colorSet(1,:),colorSet(2,:),BP,BTI);
set(gca, 'XScale', 'log', 'YScale', 'linear');
xlim(xlims_samples);
ylim(ylims_samples)
title({"Power Law with Exponential Cut-off", "\alpha=2.5, \lambda=1E-2"})
xlabel("Number of Samples")
ylabel("Burstiness Metric Value")

subplot(3,2,4)
avg = avgFun(samples(:));
disp(type + ": " + num2str(avg))
plot_windowed(samples, colorSet(1,:),colorSet(2,:),BP,BTI);




set(gca, 'XScale', 'log', 'YScale', 'linear');
xlim(xlims_window);
ylim(ylims_window)
title({"Power Law with Exponential Cut-off", "\alpha=2.5, \lambda=1E-2"})
xlabel("Time Window T_f");
ylabel("Burstiness Metric Value");


subplot(3,2,5)
type = 'lgnorm';
args = {1, 2};
mu = 1;
sigma = 2;
pd = makedist('Lognormal', 'mu', mu, 'sigma', sigma);
mean_val = mean(pd);
std_val = std(pd);
samples = random('Lognormal', mu, sigma, M_samples, N_samples);
quants = logninv(p_BTI,mu,sigma);
CoV = std_val / mean_val;
BP = (CoV - 1) / (CoV + 1);
BTI = get_tail_ratio_from_quantiles(quants,p_BTI);

plot_per_sample_size(samples, colorSet(1,:),colorSet(2,:),BP,BTI);
set(gca, 'XScale', 'log', 'YScale', 'linear');
xlim(xlims_samples);
ylim(ylims_samples)
title("Lognormal, \mu=1 \sigma=2")
xlabel("Number of Samples")
ylabel("Burstiness Metric Value")

subplot(3,2,6)
avg = avgFun(samples(:));
disp(type + ": " + num2str(avg))
plot_windowed(samples/avg, colorSet(1,:),colorSet(2,:),BP,BTI);
set(gca, 'XScale', 'log', 'YScale', 'linear');
xlim(xlims_window);
ylim(ylims_window)
title("Lognormal, \mu=1 \sigma=2")
xlabel("Time Window T_f");
ylabel("Burstiness Metric Value");

fontsize(gcf,"increase")
fontsize(gcf,"increase")
fontsize(gcf,"increase")



function plot_windowed(samples,color1,color2,BP,BTI)
    N_window = 25;
    xlims = [1E1 1E6];
    samples = samples(:);
    windowed_metrics = zeros(2,N_window);
    t_window = logspace(log10(xlims(1)),log10(xlims(2)),N_window);

    for ii = 1:N_window
        tf = linspace(0,t_window(ii),length(samples))';
        samples_windowed = min(samples,tf);
        std_emp = std(samples_windowed);
        mu_emp = mean(samples_windowed);
        windowed_metrics(1,ii) = (std_emp - mu_emp)./(std_emp + mu_emp);
        windowed_metrics(2,ii) = get_sample_BTI(samples_windowed);
    end

    semilogx(t_window,windowed_metrics(1,:),'LineWidth',3,'Color',color1);
    hold on
    semilogx(t_window,windowed_metrics(2,:),'LineWidth',3,'Color',color2);  
    semilogx(xlims,BP*ones(1,2),':','Color',color1,'LineWidth',2)
    semilogx(xlims,BTI*ones(1,2),':','Color',color2,'LineWidth',2)
    xlim(xlims)
    grid on

end
