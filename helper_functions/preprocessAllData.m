function allSegments = preprocessAllData()
% preprocessAllData  Reads all CSV files for child, deskwork, housework, rest
%   and for participants 1..5. Splits each dataset into valid segments
%   by ignoring:
%     - gaps > 60s => new segment
%     - segments with > 30 gaps > 10s => entire segment discarded
%
% Output: a struct array allSegments with fields:
%   .groupName    (string)
%   .participantData(p).validSegments, which is a cell array of segments
%       each segment has fields 't' and 'A'
%
% Example:  allSegments(1).groupName = 'child'
%           allSegments(1).participantData(1).validSegments{1} = struct('t',..., 'A',...)
%           ...
%
% EXPECTS DATA` AS FOLLOWS:
% "dataset" folder is the unzipped file from Sano and Takeuchi's public
% dataset,https://infoshako.sk.tsukuba.ac.jp/~databank/db/DB2023-02.zip
% with data description https://commons.sk.tsukuba.ac.jp/wp-content/uploads/sites/13/2023/10/dbp_2023-02en.pdf
% Takeuchi, M., & Sano, Y. (2024). “Burstiness of human physical activities
% and their characterization,” J. Comput. Soc. Sc. (2024).
% https://doi.org/10.1007/s42001-024-00247-w

    groups = {'child','deskwork','housework','rest'};
    nGroups = numel(groups);
    nParticipants = 5;

    % Pre-allocate struct array
    % We'll store for each group an array of participant data
    allSegments = struct('groupName',cell(nGroups,1), ...
                         'participantData',[]);
    
    for g = 1:nGroups
        groupName = groups{g};
        allSegments(g).groupName = groupName;
        % Each group has up to 5 participants
        participantStruct = struct('validSegments',cell(nParticipants,1));
        
        for p = 1:nParticipants
            filename = sprintf('dataset/%s%d.csv', groupName, p);

            if ~isfile(filename)
                warning('File not found: %s. Skipping...', filename);
                participantStruct(p).validSegments = {};
                continue;
            end

            % Load the data
            T = readtable(filename); % columns: t (ms), A
            time_ms = T.t;
            A = T.a;

            time_s = double(time_ms)/1000;   % convert ms->s

            % Sort by time if needed
            [time_s, sortIdx] = sort(time_s);
            A = A(sortIdx);

            % Identify big gaps > 60s => separate segments
            dt = diff(time_s);
            bigGapIdx = find(dt>60);
            disp(sum(bigGapIdx))

            segStart = [1; bigGapIdx+1];
            segEnd   = [bigGapIdx; length(time_s)];

            validSegs = {};

            for s = 1:length(segStart)
                idxRange = segStart(s):segEnd(s);
                if length(idxRange)<2
                    continue; % trivial or empty
                end

                % Check how many local gaps >10s
                local_dt = dt(idxRange(1:end-1));
                numBig10 = sum(local_dt > 10);

                if numBig10 <= 30
                    % Accept this segment
                    segT = time_s(idxRange);
                    segA = A(idxRange);
                    segStruct = struct('t', segT, 'A', segA);
                    validSegs{end+1} = segStruct; %#ok<AGROW>
                else
                    disp('too many gaps!')
                end
            end

            participantStruct(p).validSegments = validSegs;
        end % participants

        allSegments(g).participantData = participantStruct;
    end % groups
end
