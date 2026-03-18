function [childIET, deskIET, houseIET, restIET] = getIETforAlpha(allSegments, alphaVal)
   min_IET_threshold = .3;


% getIETforAlpha  Given the preprocessed segments (allSegments) and a
%   threshold alphaVal, compute the inter-event times for child, deskwork,
%   housework, and rest. We look at each segment and define E(t)=1 if
%   A(t)>alphaVal. Then compute the difference in times for consecutive 1s.
%
% Output: four arrays of inter-event times: childIET, deskIET, houseIET, restIET.

    % Initialize output
    childIET     = [];
    deskIET      = [];
    houseIET     = [];
    restIET      = [];

    for g = 1:numel(allSegments)
        groupName = allSegments(g).groupName;
        partData  = allSegments(g).participantData;

        % Collect IET for this group
        groupIET = [];

        % For each participant
        for p = 1:numel(partData)
            segs = partData(p).validSegments;
            if isempty(segs)
                continue;
            end

            thisFileIET = [];
            % For each valid segment
            for s = 1:numel(segs)
                segT = segs{s}.t;
                segA = segs{s}.A;

                E = (segA > alphaVal);

                % Indices where E(t)=1
                eventIdx = find(E==1);
                if numel(eventIdx)<2
                    continue; % not enough to form an interval
                end

                % times of the events
                eventTimes = segT(eventIdx);
                localIETs = diff(eventTimes);
                thisFileIET = [thisFileIET; localIETs]; %#ok<AGROW>
            end

            thisFileIET = thisFileIET(thisFileIET>=min_IET_threshold);
            groupIET = [groupIET; thisFileIET]; %#ok<AGROW>
        end

        % Assign to the correct output array
        switch groupName
            case 'child'
                childIET = groupIET;
            case 'deskwork'
                deskIET  = groupIET;
            case 'housework'
                houseIET = groupIET;
            case 'rest'
                restIET  = groupIET;
        end
    end
end
