% function that gets the indices of the fUS block (.bin) files (stored in
% the fUS timelog.txt file) that correspond to the attempted trials in the
% behavioral data. 
% Uses session info, behavior log, and session/run combination. 
% returns the fUS block indices, the modified behavior file (the trials are
% cut off to reflect only the timing where fUS acquisition was running),
% and the trial Length in seconds. Note that trial length includes a buffer
% period past the end of the trial so we can see evolution of hemodynamics.
% This is not reflected in the bare behavior file, and that is why we pass
% it back here.
%
% Further details of function output
% trialLength - int; [s]
% 
% [fUSBlockIndex, behavior, trialLength] = getfUSLogIndices(SessionInfo, behavior, Session, Run)

function [fUSBlockIndex, behavior, trialLength] = getfUSLogIndices(SessionInfo, behavior, Session, Run)
if nargin < 4
    error('this function requires all arguments. type ''help getfUSLogIndices'' for function definition')
end

% get the row of interest so we can check if this run had labjack support
sessionRows = find(SessionInfo.Session == Session);
runRows = find(SessionInfo.Run == Run);
rowOfInterest = intersect(sessionRows,runRows);

%% find trial start based on memory start times 

% get fixation, memory, and movement epochs
bufferPeriod = 5;   % [s] (buffer period into the ITI to see time traces more clearly)
[fix, mem, mov] = getEpochs(behavior); % this is where we actually get the epoch lengths!!!
fixationPeriodLength = ceil(fix(2)-fix(1));
trialLength = ceil(mov(end)-fix(1)+bufferPeriod); % Round up to closest int

fprintf('\nFixation period set to %i sec before memory delay period onset\n',fixationPeriodLength)
fprintf('ITI buffer period set to %i sec after reward is dispensed\n',bufferPeriod)
fprintf('Total trial length = %2.1f\n\n',trialLength)

% get trial times from behavior, aligning to memory start times
if SessionInfo.Labjack(rowOfInterest) && (Session<30 || Session == 35 || Session == 39)   
    t_memstarts = ...
        labjackCorrection(Session,Run,behavior);    %find memory start times (labjack corrected)
else
    % for newer time logs, we don't want to use labjack correction (even if
    % its supported & data exists) because we will pull the tcp/ip sync'd
    % data below (where fUS time logs are based on behavior pc clock).
    try
        t_memstarts = parseTime(behavior,'memory');     %find memory start times (behavior data only)
    catch
        warning('No memory period found. Using cue start instead');
        t_memstarts = parseTime(behavior,'cue');
end
t_trialstarts = t_memstarts-fixationPeriodLength;         %calculate trial start times

%% find fUS acqs start times (Verasonics' log file)
cd(getUserDataPath('pathType','fUS','session',Session,'run',Run));
t_fUSacqs = abs(textread('timeLog.txt', '%f')); %#ok<DTXTRD> % load fUS acq time log

% reshape fUS time log to [HH MM SS.MS]
if Session < 30 || Session == 35 || (Session == 39 && Run ~= 1)
    % older style fUS time log has 4 columns
    t_fUSacqs = reshape(t_fUSacqs,4, [])';  
    t_fUSacqs = t_fUSacqs(:,1:3);                   
else
    % new fUS time log has 7 columns, 4-6 are already sync'd to behavior ctrlr
    t_fUSacqs = reshape(t_fUSacqs,7, [])';
    t_fUSacqs = t_fUSacqs(:,4:6);                   
end

T_fUS = t_fUSacqs(:,1)*3600 + t_fUSacqs(:,2)*60+ t_fUSacqs(:,3); % convert to [s]

%% Truncate behavior time log if t_starttrials(i) > t_fUS(end)-epochLength

% finds index of last trial that was mapped by fUS
iCut = find((t_trialstarts+trialLength) < T_fUS(end), 1, 'last');

% Truncate behavior if it kept going after fUS recording
if ~isempty(iCut)
    behavior(iCut:end) = [];
    t_trialstarts(iCut:end) = [];
end

% finds index of first trial that was mapped by fUS
iCut = find((t_trialstarts+fix(1)) < T_fUS(1), 1, 'last');
% truncate behavior if it started before fUS recording
if ~isempty(iCut)    
    behavior(1:iCut) = [];
    t_trialstarts(1:iCut) = [];
end    

% just to be double-sure
if (t_trialstarts(1) < T_fUS(1)) || (t_trialstarts(end) > T_fUS(end))
    error('behavior data exists outside fUS recording times!')
end

%% matching time logs (find first t_fUSacq before t_trialstart)
nTrials = size(behavior,2);
% fUS indices for left cue trials
fUSBlockIndex = NaN(nTrials,1);
for i = 1:nTrials    
    [~, fUSBlockIndex(i)] = min(abs(T_fUS - t_trialstarts(i)));
end

%% removing errant values 
% this is necessary because, in some cases, the behavior may have kept
% running while the fUS acquisition was stopped & re-started (see S7R1). 
[fUSBlockIndex, uniqueIndices] = unique(fUSBlockIndex);
behavior = behavior(uniqueIndices);

end