function [fix, mem, mov] = getEpochs(behavior)
% [fix, mem, mov] = getEpochs(behavior)
% returns the epoch/phase timing of the trial based on the behavior
% example:
% behavior is a 1xnTrials struct with fields, e.g. trialID, trialstart,...
% fix = [-5.2 0]
% mem = [0  4.9]
% mov = [4.9 12]
% note: all values are in seconds relative to initial/memory cue onset

% make sure we didn't send TBI data to this NHP function by accident.
if istable(behavior)
    warning('this is TBI behavior data. Calling the appropriate function.')
    [fix, mem, mov] = TBI.getEpochs(behavior); 
    return
elseif isfield(behavior,'RewardTime')
    warning('this is TBI behavior data. Calling the appropriate function.')
    [fix, mem, mov] = TBI.getEpochs(behavior);
    return
end

% only look at successful trials
success = cell2mat({behavior.success});
behavior(~success) = [];

% get vector of times for each phase
fixationHold = parseTime(behavior, 'fixationhold');
cue = parseTime(behavior, 'cue');
try
    memoryStart = parseTime(behavior, 'memory');
    AlignToTime = 'memory';
catch
   warning('No memory period. Aligning to cue period instead.');
   AlignToTime = 'cue';
end
        
goCue = parseTime(behavior,'target_acquire');
reward = parseTime(behavior, 'reward');

% realign to memory or cue start
switch AlignToTime
    case 'memory'
        fixationHold = fixationHold - memoryStart;
        cue = cue - memoryStart;
        goCue = goCue - memoryStart;
        reward = reward - memoryStart;
    case 'cue'
        fixationHold = fixationHold - cue;
        goCue = goCue - cue;
        reward = reward - cue;
end

% set periods
fix = [median(fixationHold) 0];
mem = [0 median(goCue)];
mov = [median(goCue) max(reward)];


%% secondary function to get vector of epoch times
function numTime = parseTime(chartime, columnTitle) 
    chartime = eval(['cell2mat({chartime.' columnTitle '}'');']);

    % extract hours/min/sec/ms from chartimes
    %#ok<*ST2NM>
    HH = str2num(chartime(:,1:2));   % Hours array
    MM = str2num(chartime(:,4:5));   % Minutes array
    SS = str2num(chartime(:,7:8));   % seconds array
    MS = str2num(chartime(:,10:12)); % milliseconds array

    % concat visual stim trial times [HH MM SS.MS]
    numTimeArray = cat(2, HH, MM, SS+1e-3*MS);

    % convert time to [s]
    numTime =   numTimeArray(:,1)*3600 + numTimeArray(:,2)*60+ numTimeArray(:,3);
end
end