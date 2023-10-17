function [t, parameterTable, parameter_counter, cnt_attempt] = MemoryMovementTask(CurrentSession, ...
    parameterTable, ...
    parameter_counter, ...
    cnt_attempt, ...
    session, ...
    overlap)
%% Parses behavioral data for the memory & overlap saccade tasks
%
% Inputs: 
% CurrentSession: Contains all the behavioral data to be parsed
% parameterTable: Proper dimensioned table for the parameters
% parameter_counter: Current record of parameters across entire datafile
% cnt_attempt: Current record of attempted trials across entire datafile
% session: session number for the waitbar.
% overlap: OPTIONAL boolean (turn true if this is an overlap task)
%
% Outputs:
% t: Consolidated table of behavioral information about each trial
% parameterTable: Returns record of all parameter values
% parameter_counter: Returns record of parameters across entire datafile
% cnt_attempt: Returns record of attempted trials across entire datafile.
%
% Written by WS Griggs & SL Norman

%% handle varargin (should replace with input parser if this grows)
if ~exist('overlap','var')
    overlap = false;
end

%% Function content

%find the beginning of trials
trial_startLine = find(strcmp(CurrentSession(:,2),'trial'));

%Initialization
wb = waitbar(0,sprintf('Parsing Behavior: Session: %d...',session));

% Number of trials in session
nTrials = length(trial_startLine);

% Iterate over each trial
for trial = 1:nTrials
    % Advance the waitbar once notch per trial
    waitbar(trial/nTrials,wb,sprintf('Parsing Behavior: Session: %d...',session));
    
    % data from file, but only for the current trial
    if trial ~= nTrials
        CurrentTrial = CurrentSession( trial_startLine(trial): trial_startLine(trial+1)-1,:);
    else
        CurrentTrial = CurrentSession(trial_startLine(trial):end,:);
    end
    
    % Information about states in current trial
    state_startLine = find(strcmp(CurrentTrial(:,2),'state'));
    nStates = length(state_startLine);
    
    % Iterate over each state
    for state = 1:nStates
        
        % Get data from current state
        if state ~= nStates
            CurrentState = CurrentTrial( state_startLine(state): state_startLine(state+1)-1,:);
        else
            CurrentState = CurrentTrial(state_startLine(state):end,:);
        end
        
        currentState_cell = CurrentState(strcmp('targetname',CurrentState(:,2)),3);
        
        % If things happened during the state, then parse it
        if ~isempty(currentState_cell)
            % Clear states_struct to prevent information from previous
            % states being carried over.
            clear states_struct;
            
            % Convert cell array into struct array
            for i = 1:length(currentState_cell)
                states_struct(i,1) = currentState_cell{i,1};
            end
            
            % iterate over each target within each state
            targetNames = unique({states_struct.targetname});
            for target = 1:length(targetNames)
                % Iterate over each parameter within each target
                parameterNames = {states_struct(strcmp({states_struct.targetname},targetNames(target))).parameter};
                for parameter = 1:length(parameterNames)
                    % Fetch session number
                    parameterTable.sessionNum(parameter_counter) = session;
                    % Fetch trial number
                    parameterTable.trialNum(parameter_counter) = str2num(CurrentSession{trial_startLine(trial),3});
                    % Fetch state name
                    parameterTable.state(parameter_counter) = CurrentTrial(state_startLine(state),3);
                    % Fetch target name
                    parameterTable.target(parameter_counter) = targetNames(target);
                    % Fetch parameter names
                    parameterTable.parameter(parameter_counter) = parameterNames(parameter);
                    % Fetch value of the parameter
                    parameterTable.value(parameter_counter) = {states_struct(strcmp({states_struct.targetname},targetNames(target)) ...
                        & strcmp({states_struct.parameter},parameterNames(parameter))).value};
                    % Advance parameter counter by one
                    parameter_counter = parameter_counter+1;
                end % End parameter loop
            end % End target loop
        end % End if statement
    end % End state loop
    
    
    % Get information about targets
    targets_cell = CurrentTrial(strcmp('target',CurrentTrial(:,2)),3);
    % Clear targets struct to prevent overlap across states
    clear targets_struct
    for i = 1:length(targets_cell)
        targets_struct(i,1) = targets_cell{i};
    end
    
    % Get information about images (if any)
    targets_cell = struct2cell(targets_struct)'; % Removes structs from earlier targets_cell
    fieldnames = fields(targets_struct);
    
    for target = 1:size(targets_cell,1)
        targetname = targets_cell{target,strcmp(fieldnames,'targetName')};
        targets_struct(target,1).image = getParameterValue(parameterTable,...
            {session, trial, 'cue',targetname,'image_number'});
    end
    
    
    % set the state names of interest in the task
    if overlap
        states = {'initialfixation', 'fixationhold', 'cue',  ...
            'target_acquire', 'target_hold', 'reward'};
        attempt_cutoff = 4; 
        success_cutoff = 6;
    else
        states = {'initialfixation', 'fixationhold', 'cue', 'memory'  ...
            'target_acquire', 'target_hold', 'reward'};    
        attempt_cutoff = 5;
        success_cutoff = 7;
    end
    
    % attempted (monkey got to cue) (>= fix, hold, cue, [memory], iti, ...)
    if nStates >= attempt_cutoff 
        
        cnt_attempt                = cnt_attempt + 1;
        
        % Extract trial number
        t(cnt_attempt).trialID     = str2num(CurrentSession{trial_startLine(trial),3});
        
        % set t(trial) phase/state entry times
        t(cnt_attempt).trialstart        = TimeTrialStarts(CurrentTrial{strcmp('trial',CurrentTrial(:,2)),1});
        for state = 1:nStates-1
            eval(['t(cnt_attempt).' states{state} ' = TimeTrialStarts(CurrentTrial{state_startLine(state),1});'])
        end
        % skip to the iti
        t(cnt_attempt).iti = TimeTrialStarts(CurrentTrial{state_startLine(nStates),1});
        
        % set other attributes
        if size(targets_struct,1)>1
            t(cnt_attempt).free_choice       = ~isequal(targets_struct(1).targetPos,targets_struct(2).targetPos);
        else
            t(cnt_attempt).free_choice = false;
        end
        t(cnt_attempt).target            = sign(targets_struct(1).targetPos(1)); %Whether target is on left or right of screen midline
        t(cnt_attempt).targetPos         = targets_struct(1).targetPos; % Actual target position
        t(cnt_attempt).success           = nStates > success_cutoff; % Successful trial if it makes to at least n states (fix, hold, cue, [memory], targetAcq, hold, reward, iti)
        t(cnt_attempt).lastState         = states{nStates-1}; % Minus one because we don't count ITI state
        t(cnt_attempt).effector          = CurrentTrial{strcmp('effector',CurrentTrial(:,2)),3}.effector; % Effector used in the task
        t(cnt_attempt).constrainedEffector = CurrentTrial{strcmp('effector',CurrentTrial(:,2)),3}.windowconstraint; %Effector that was restricted in the task
        t(cnt_attempt).targetInfo        = targets_struct; % Total information about target. Partially redundant with targetPos variable
        
        % If free choice and successful trial, indicate which target was
        % chosen
        if t(cnt_attempt).free_choice && t(cnt_attempt).success && any(strcmp('choice',CurrentTrial(:,2)))
           ImageChosen = CurrentTrial{strcmp('choice',CurrentTrial(:,2)),3};   
           ImageNum = str2num(ImageChosen(length('image ')+1:end));
           t(cnt_attempt).chosenTarget =  targets_struct([targets_struct.image] == ImageNum,:);
        end
        
        % If there a reported RWD size, save it. Otherwise, use NaN.
        RWD_ind = find(strcmp(CurrentTrial(:,2),'RWD magnitude'));
        if RWD_ind
            t(cnt_attempt).RWDsize       = str2num(cell2mat(CurrentTrial(RWD_ind,3)));
        else
            t(cnt_attempt).RWDsize       = NaN;
        end
        
        % Output the actual reward delivered on successful trials.
        % This should agree with RWDsize, but it may differ slightly.
        if t(cnt_attempt).success
           ITI_start = datenum(t(cnt_attempt).iti,'HH:MM:SS,FFF');
           reward_start = datenum(t(cnt_attempt).reward,'HH:MM:SS,FFF');
           t(cnt_attempt).RWD_delivered = ITI_start - reward_start;
           
           % Convert from fraction of a day into ms
           t(cnt_attempt).RWD_delivered = t(cnt_attempt).RWD_delivered * 24 * 60 * 60 * 1000;
        end
        
        % If there is a reported blocknum, save it. Otherwise, use NaN
        blocknum_ind = find(strcmp(CurrentTrial(:,2),'blocknum'));
        if blocknum_ind
            t(cnt_attempt).blocknum     = str2num(cell2mat(CurrentTrial(blocknum_ind,3)));
        else
            t(cnt_attempt).blocknum     = NaN;
        end
    end %End nStates loop
end % End trial loop
close(wb);

% re order field for easier legibility
if overlap
    t = moveField(t, 'cue', 'fixationhold');
    t = moveField(t, 'target_acquire', 'cue');
else
    t = moveField(t, 'target_acquire', 'memory');
end
t = moveField(t, 'target_hold', 'target_acquire');
t = moveField(t, 'reward', 'target_hold');
end

function Start_Time = TimeTrialStarts(TimeString)
%%
% Written by Sumner (I think)
id_time_start     = strfind(TimeString,' ');
Start_Time        = TimeString(id_time_start + 1:end);
end

