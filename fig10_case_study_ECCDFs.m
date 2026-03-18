figure(992)
clf

alphaVals = [1 15 100];
colors = lines(length(alphaVals));
allSegments = preprocessAllData();

data_labels = {'Child Play','Desk Work','Housework','Rest'};
panel_labels = {'(A)','(B)','(C)','(D)'};
N_datasets = length(data_labels);

for ii = 1:length(alphaVals)

    alph = alphaVals(ii);
    % childIET, deskIET, houseIET, restIET
    [data{1:N_datasets}] = getIETforAlpha(allSegments, alph);

    for jj = 1:N_datasets
        subplot(2,2,jj);
        [f,x] = ecdf(data{jj});

        mu = mean(data{jj});
        %minx = min(data{jj});
        x_scaled = x/mu;
        hold on
        loglog(x_scaled, 1-f, 'Color', colors(ii,:),'LineWidth',3,'DisplayName', "\theta = " +  alph + " m/s^3");
    end
end

%sgtitle("ECCDFs of IETs of Human Activity")
fontsize("increase")
fontsize("increase")
fontsize("increase")

for jj = 1:N_datasets
    subplot(2,2,jj);   
    hold on
    loglog(x_scaled, expcdf(x_scaled,"upper"),'g:','LineWidth',3,'DisplayName',"exponential")
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    grid on

    ylabel('P[X]>x')
    xlabel('IET (normalized)')
    %xlim([min(data{jj}) max(data{jj})])
    %mu = mean(data{jj});
    %xlim([1 max(data{jj})/mu])
    xlim([0.5 1E2])
    ylim([1E-4 1])
    title(panel_labels{jj})
    legend()
end

fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")