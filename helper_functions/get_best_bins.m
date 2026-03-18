function [binEdges,logBinCounts] = get_best_bins(counts,isData)
    % Written by ChatGPT 4.0 on 02/01/2024


    if ~exist('isData','var')
        isData = false;
    end

    if isData
        data = counts;
        min_val = min(data) - 0.001;
        max_val = max(data) + 0.001;
    else
        min_val = find(counts,1,"first");
        max_val = find(counts,1,"last");
    end



    if ~exist('startBinNum','var')
        startBinNum = 50;
    end
    % Initialize logarithmic bin edges
    binEdges = logspace(log10(min_val), log10(max_val), startBinNum); % Adjust the number of bins as needed
    %binEdges = [binEdges, inf]; % Append infinity to cover the upper edge of the last bin
    
    % Initialize the array to hold the new bin counts
    logBinCounts = zeros(1, length(binEdges)-1);
    
    % Perform iterative re-binning
    while true

        if isData
            [logBinCounts,~] = histcounts(data,binEdges);
        else
            for ii = 1:length(counts)
                % Find the bin index for the current value
                binIndex = find(ii < binEdges, 1) - 1;
                %disp(binIndex)
                if isempty(binIndex)
                    binIndex = length(logBinCounts);
                end
                % Add the count to the appropriate bin
                logBinCounts(binIndex) = logBinCounts(binIndex) + counts(ii);
            end
        end
        
        % Check for empty bins (after the first few)
        if any(logBinCounts(1:end) == 0)
            % Reduce the number of bins
            binEdges = logspace(log10(min_val), log10(max_val), length(binEdges)-1);

            if length(binEdges) <= 5
                disp("can't do it!")
                return
            end

            %binEdges = [binEdges, inf]; % Append infinity to cover the upper edge of the last bin
            logBinCounts = zeros(1, length(binEdges)-1);
        else
            break; % Stop if there are no empty bins
        end
    end

end