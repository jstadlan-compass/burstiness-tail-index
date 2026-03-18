
% Plot BTI vs. BP for the activity thresholds 
figure(332)
clf
tiledlayout(1,2,'TileSpacing','Compact');
nexttile
scatterplot_BTI_BP_ranges(0:15, 1)
title("\alpha = 1 - 15 m/s^3")
a = colorbar;
a.Label.String = '\alpha (m/s^3)';
xlabel("BP")
ylabel("BTI")
nexttile
scatterplot_BTI_BP_ranges(30:20:360, 1)
a = colorbar;
a.Label.String = '\alpha (m/s^3)';
title("\alpha = 30 - 350 m/s^3")
xlabel("BP")
ylabel("BTI")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")
fontsize("increase")

sgtitle("Human Activity IETs at Varying Jerk Thresholds")


function [resultTableCV, resultTableBTI] = scatterplot_BTI_BP_ranges(alphaVec, chunkSize)
    
    addpath("burstiness_helper_functions")

    % 1) Load/parse all data once
    allSegments = preprocessAllData();

    % 2) Choose the alpha thresholds
    % small
    %alphaVec = 1:15;
    %chunkSize = 1;

    % large
    %alphaVec = 1:300;
    %chunkSize = 25;

    nAlpha = length(alphaVec);

    % Initialize matrices for CV and BTI
    CV_results = nan(nAlpha,4);   % columns: child, deskwork, housework, rest
    BTI_results = nan(nAlpha,4); % columns: child, deskwork, housework, rest

    % 3) Loop over alpha values, compute CV and BTI for each group
    for iAlpha = 1:nAlpha
        alphaVal = alphaVec(iAlpha);

        [childIET, deskIET, houseIET, restIET] = getIETforAlpha(allSegments, alphaVal);

        CV_child     = computeCV(childIET);
        CV_deskwork  = computeCV(deskIET);
        CV_housework = computeCV(houseIET);
        CV_rest      = computeCV(restIET);

        BTI_child     = get_tail_ratio(childIET);
        BTI_deskwork  = get_tail_ratio(deskIET);
        BTI_housework = get_tail_ratio(houseIET);
        BTI_rest      = get_tail_ratio(restIET);

        CV_results(iAlpha,:)  = [CV_child,  CV_deskwork,  CV_housework,  CV_rest];
        BTI_results(iAlpha,:) = [BTI_child, BTI_deskwork, BTI_housework, BTI_rest];
    end

    % 4) Create and display tables for CV and BTI
    resultTableCV = table(alphaVec(:), ...
        CV_results(:,1), CV_results(:,2), CV_results(:,3), CV_results(:,4), ...
        'VariableNames', {'alpha','CV_child','CV_deskwork','CV_housework','CV_rest'});
    disp(resultTableCV);

    resultTableBTI = table(alphaVec(:), ...
        BTI_results(:,1), BTI_results(:,2), BTI_results(:,3), BTI_results(:,4), ...
        'VariableNames', {'alpha','BTI_child','BTI_deskwork','BTI_housework','BTI_rest'});
    disp(resultTableBTI);
 %% 4) Group alpha into non-overlapping ranges of length 5

    % In this example: 1–5, 6–10, 11–15, 16–20, 21–25, 26–30
    nChunks = nAlpha / chunkSize;  % = 30/5 = 6

    chunkedCV  = nan(nChunks,4);
    chunkedBTI = nan(nChunks,4);
    chunkAlpha = nan(nChunks,1);   % We'll store the midpoint alpha for color

    for c = 1:nChunks
        % Indices in alphaVec for this chunk
        idxStart = (c-1)*chunkSize + 1;
        idxEnd   = c*chunkSize;
        theseIdx = idxStart:idxEnd;

        % Average across each group’s values in that alpha range
        chunkedCV(c, :)  = mean(CV_results(theseIdx, :), 1, 'omitnan');
        chunkedBTI(c, :) = mean(BTI_results(theseIdx, :),1, 'omitnan');

        % We can color by the *midpoint* alpha of this chunk
        chunkAlpha(c) = mean(alphaVec(theseIdx));
    end

    %% 5) Plot (average CV vs. average BTI), colored by chunkAlpha, shaped by group
    hold on;
    markerShapes = {'o','s','^','d'};  % child, deskwork, housework, rest
    groupLabels  = {'child','deskwork','housework','rest'};

    colormap winter

    for g = 1:4
        % X-values = chunkedCV(:, g), Y-values = chunkedBTI(:, g)
        scatter(...
            chunkedCV(:,g), ...
            chunkedBTI(:,g), ...
            400, ...                % marker size
            chunkAlpha, ...        % color data
            markerShapes{g}, ...   % shape depends on group
            'filled' ...
        );
    end

    colorbar;   
    clim([min(chunkAlpha), max(chunkAlpha)]);
    xlabel('BP (moving average)');
    ylabel('BTI (moving average)');
    legend(groupLabels, 'Location','best');
    title('BTI vs. BP over varying thresholds');
    hold off;
    xlim([-0.5 1])
    ylim([-0.5 1])
    axis square
    grid on


end

% ========================================================================

% ========================================================================

% ========================================================================
function CV = computeCV(x)
% computeCV  Returns std(x)/mean(x), or NaN if x is empty/degenerate.
    if isempty(x) || all(x == x(1))
        CV = NaN;
        return;
    end
    CV = burstiness(x);
end

