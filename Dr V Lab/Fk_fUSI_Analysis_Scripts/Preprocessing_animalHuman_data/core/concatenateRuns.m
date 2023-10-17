function [dopOut, LIndsOut, RIndsOut, behaviorOut, coreParams, wallpaper, UFout]...
    = concatenateRuns(sessionRunList, registered)
% [dopOut, LIndsOut, RIndsOut, behaviorOut, coreParams] = concatenateRuns(sessionRunList)
%
% session/run list is of the format:
% session run
% [  1     1  ]
% [  1     2  ]
% [  2     1  ]
% [    ...    ]
%
% registered is a boolean whether you want to request motion corrected data

% When concatenating registered data, this function does not currently
% align between sessions, i.e. there is quite likely a position shift
% between sessions.

%% check that data type is consistent across the runs requested
load('SessionInfo.mat','SessionInfo');

% Preallocate space
indx = NaN(size(sessionRunList,1),1);

for i = 1:size(sessionRunList,1)
    % set current session/run combo & find row in the SessionInfo table
    session = sessionRunList(i,1);
    run = sessionRunList(i,2);
    indx(i) = intersect(find(SessionInfo.Session==session),...
        find(SessionInfo.Run==run));
end

% Check if more than one monkey is used
if length(unique(SessionInfo.Monkey(indx)))>1
    monkeyPass = false;
else
    monkeyPass = true;
end

% Check if all are same task
if length(unique(SessionInfo.Task(indx)))>1
    taskPass = false;
else
    taskPass = true;
end

% Check if same number of probes was used consistently
if any(strcmp(SessionInfo.Probe(indx),'DoubleProbe'))
    if any(~strcmp(SessionInfo.Probe(indx),'DoubleProbe'))
        probePass = false;
    else
        probePass = true;
    end
else
    probePass = true;
end

% Pause if the sessions are deemed incompatible
if ~(monkeyPass && taskPass && probePass)
    warning('The selected sessions are not compatible!')
    input('Press any key to continue. (ctrl+c to quit)')
end

%% load the data and concatenate
for i = 1:size(sessionRunList,1)
    % create filename string for this session/run
    session = sessionRunList(i,1);
    run = sessionRunList(i,2);
    dopplerFileName = strcat('doppler_S', num2str(session), '_R', num2str(run), '.mat');
    
    % Only load motion corrected data if not concatenated.
    if registered && size(sessionRunList,1)==1
        dopplerFileName = strrep(dopplerFileName,'.mat','+normcorre.mat');
    end
    
    % load data
    try
        try
            fprintf('Loading ''%s'' from %s...\n',dopplerFileName,getUserDataPath('pathType','Processed','speed','fast'))
            load(fullfile(getUserDataPath('pathType','Processed','speed','fast'), dopplerFileName), 'iDop', 'LInds', 'RInds', 'behavior', 'coreParams','UF','wallpaper')
        catch
            warning('Data not found in "fast" directory. Looking elsewhere.');
            fprintf('Loading ''%s'' from %s...\n',dopplerFileName,getUserDataPath('pathType','Processed'))
            load(fullfile(getUserDataPath('pathType','Processed'), dopplerFileName), 'iDop', 'LInds', 'RInds', 'behavior', 'coreParams','UF','wallpaper')
        end
    catch
        if registered
            warning('Motion corrected data not found. Attempting to load unregistered data.')
            dopplerFileName = strrep(dopplerFileName,'+normcorre','');
            try
                fprintf('Loading ''%s'' from %s...\n',dopplerFileName,getUserDataPath('pathType','Processed','speed','fast'))
                load(fullfile(getUserDataPath('pathType','Processed','speed','fast'),dopplerFileName), 'iDop', 'LInds', 'RInds', 'behavior', 'coreParams','UF');
            catch
                warning('Unregisterd data not found in "fast" directory. Looking elsewhere.');
                fprintf('Loading ''%s'' from %s...\n',dopplerFileName,getUserDataPath('pathType','Processed'))
                load(fullfile(getUserDataPath('pathType','Processed'), dopplerFileName), 'iDop', 'LInds', 'RInds', 'behavior', 'coreParams','UF')
            end
        else
            warning('Cannot load specified file.');
        end
    end
    
    if i==1
        % start by initializing to first instance
        dopOut = iDop;
        LIndsOut = LInds;
        RIndsOut = RInds;
        behaviorOut = behavior;
        UFout = UF;
    else % we are on the second run or later, so concatenate them
        % make sure that we have the same trial length (sometimes, if a run
        % has an epoch length that jitters above or below a full second,
        % this can get rounded up/down based on jitter between runs.
        
        if ~all([size(dopOut,1)==size(iDop,1),size(dopOut,2)==size(iDop,2),size(dopOut,3)==size(iDop,3)])
            minDims = min(size(dopOut), size(iDop));
            dopOut = dopOut(1:minDims(1), 1:minDims(2), 1:minDims(3), :);
            iDop = iDop(1:minDims(1), 1:minDims(2), 1:minDims(3), :);
        end
        % concatenate across trials (4th dim)
        tmp = cat(4, dopOut, iDop);
        dopOut = tmp;
        % concatenate L/R indices
        LIndsOut = cat(2, LIndsOut, LInds);
        RIndsOut = cat(2, RIndsOut, RInds);
        
        % update UF       
        if ~isfield(UF,'depth')
            warning('UF missing depth field')
        elseif ~isequal(UF.Depth,UFout.Depth)
            warning('Different depths between the sessions. Using depth from first run');
        end
        % update behavior
        try
            behaviorOut = [behaviorOut behavior];
            fprintf('Successful trials concatenated from %i to %i\n', ...
                size(behaviorOut,2), size(behaviorOut,2)+size(behavior,2))
        catch
            disp('Could not concatenate behavior files');
        end % Try catch loop for behavior
    end % If over session number
end % For loop over runs

% save out a pre-image-registration wallpaper
if ~exist('wallpaper','var')
    wallpaper = nthroot(squeeze(mean(mean(dopOut,3),4)),3);
    lowerCutoff = quantile(wallpaper(:),0.01);
    wallpaper(wallpaper<lowerCutoff) = lowerCutoff;
end

% Motion correct entire concatenated dataset. Needs to be done at
% the very end to account for shifts in position between sessions.
if registered && size(sessionRunList,1)>1
    coreParams.motionCorrection = true;
    coreParams.method = 'normcorre';
    [dopOut, coreParams] = correctMotion(dopOut, coreParams);
end

fprintf('Passing back concatenated doppler data\n')


end %Main function end