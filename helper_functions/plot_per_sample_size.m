function plot_per_sample_size(samples,color1,color2,BP,BTI,CALC_MEAN)
    QLO = .25/2;
    QHI = 1-QLO;

    if ~exist('CALC_MEAN','var')
        CALC_MEAN = false;
    end

    Nsamples = size(samples,2);
    sample_subset_n = ceil(logspace(log10(5),log10(Nsamples),50));
    BTI_subset = zeros(size(samples,1),length(sample_subset_n));
     % get means:
    mu_samples = cumsum(samples,2)./(1:Nsamples);
    stddev_samples = sqrt(cumsum(samples.^2,2)./(1:Nsamples)-mu_samples.^2);
    CoV_samples = stddev_samples./mu_samples;

    BP_samples = (CoV_samples-1)./(CoV_samples+1);

    % Kim & Jo finite correction:
    n = repmat(1:Nsamples,size(samples,1),1);
    r = CoV_samples;
    
    BP_samples = (sqrt(n+1).*r-sqrt(n-1))./((sqrt(n+1)-2).*r+sqrt(n-1));
    BP_samples(:,1)=BP_samples(:,2);

    % Error (not normalizing because -1 to 1 bounded scale)
    %BP_samples = BP_samples-BP;

 
    
    %BP_mean = mean(BP_samples,1);
    %plot(BP_mean,'LineWidth',3)
    hold on
    %plot(BP_samples')
    
    x = 1:(Nsamples);

    if CALC_MEAN
        plot(x,mean(BP_samples,1)','LineWidth',3,'Color',color1);
    else
        plot(x,median(BP_samples,1)','LineWidth',3,'Color',color1);
    end
    hold on
    y25 = quantile(BP_samples,QLO,1);
    y75 = quantile(BP_samples,QHI,1);
    patch([x flip(x)], [y25 flip(y75)],color1, 'FaceAlpha',0.3, 'EdgeColor','none')
    grid on
    plot([1 Nsamples],BP*ones(1,2),':','Color',color1,'LineWidth',2)



    for ii = 1:size(samples,1)
        for jj = 1:length(sample_subset_n)
            N = sample_subset_n(jj);
            BTI_subset(ii,jj) = get_tail_ratio(samples(ii,1:N));
        end
    end

    % Set to 0
    %BTI_subset = BTI_subset-BTI;

    x = sample_subset_n;
    if CALC_MEAN
        plot(x,mean(BTI_subset,1),'LineWidth',3,'Color',color2);
    else
        plot(x,median(BTI_subset,1),'LineWidth',3,'Color',color2);
    end

    y25 = quantile(BTI_subset,QLO,1);
    y75 = quantile(BTI_subset,QHI,1);
    patch([x flip(x)], [y25 flip(y75)],color2, 'FaceAlpha',0.3, 'EdgeColor','none')
    plot([1 Nsamples],BTI*ones(1,2),':','Color',color2,'LineWidth',2)
    grid on
end