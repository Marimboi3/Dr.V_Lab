% This script parses the behavioral data from the .dat files into a usable
% trial marker log, t_run_#. It stores the pathing and run data in
% SessionDataPaths.mat, a small database to track the data we accumulate
% across multiple monkeys, tasks, people, and their computers.
%
% optional vargin: missingRunRows is a vector of missing rows in
% SessionInfo that can be manually defined. Otherwise, this function will
% search for missing parsed behavioral data automatically.
%
%States (successful trials):
%1. init_fixation_acquire
%2. init_fixation_hold
%3. cue
%4. memory
%5. target_acquire
%6. target_hold
%7. reward
%8. success_iti
%
% Initially written by Vasileios. Modified by Whitney and Sumner

function ParseBehavior(missingRunRows)

% load session info
tmp = load('SessionInfo.mat'); SessionInfo = tmp.SessionInfo; clear tmp
nEntries = size(SessionInfo,1);
fprintf('SessionInfo.mat loaded successfully. %i entries found.\n\n', nEntries)
startPath = pwd;

% run through the entire SessionInfo file to look for missing data
if ~exist('missingRunRows','var')
    missingRunRows = [];
    for row = 1:size(SessionInfo,1)
        
        cd(getUserDataPath('pathType','Behavioral','session',SessionInfo.Session(row),'run',SessionInfo.Run(row)));
        
        % make sure the data paths are not on path (messes with exist func)
        warning('off'), rmpath(genpath('.')), warning('on')
        
        % if we can't find the parsed behavior file where it belongs...
        if ~exist(SessionInfo.behaviorFileName(row),'file')
            missingRunRows = [missingRunRows row];
            fprintf('Session %i Run %i parsed behavior log missing.\n',...
                SessionInfo.Session(row), SessionInfo.Run(row))
        else
            fprintf('Session %i Run %i parsed behavior found!\n',...
                SessionInfo.Session(row), SessionInfo.Run(row))
        end
    end
    cd(startPath)
    
    fprintf('%i out of %i parsed behavior files found!\n\n\n',...
        nEntries - length(missingRunRows), nEntries)
    if ~isempty(missingRunRows)
        uiwait(msgbox('Click OK to generate missing files (ctrl+c to quit)'),120);
    end
else
    fprintf('Missing run row = ')
    disp(num2str(missingRunRows))
end

%% run through any/all rows that are missing parsed data
for row = missingRunRows
    
    Session = SessionInfo.Session(row);
    Run = SessionInfo.Run(row);
    fprintf('Now attempting to parse data for Session %i Run %i\n', Session, Run);
    
    % Load behavioral Data from .dat file
    cd(getUserDataPath('pathType','Behavioral','session', Session,'run', Run));
    
    Date = SessionInfo.Date(row);
    
    disp('Reading in file. Please wait...')
    if datenum(Date, 'dd-mmm-yyyy') >= datenum('22-Jan-2020','dd-mmm-yyyy')
        disp('Using new behavioral parser. ')
        rawfile = fopen(char(SessionInfo.behavioralRawFile(row))); % Open file
        
        % Initial parse of file into appropriate columns (time and message)
        file = textscan(rawfile, '%q %*q %q', ...
            'Delimiter',{' - '},'EndOfLine','\n');
        fclose(rawfile);
        file = horzcat(file{:});
        % now run the parse sub-function
        behavior = parseDatafile(file);  % (entire workspace in this scope gets saved)
    else
        % Use old data log parser
        disp('Using old behavioral parser. ')
        file = textread(char(SessionInfo.behavioralRawFile(row)), ...
            '%s', 'delimiter', '\n', 'whitespace', ''); %#ok<DTXTRD>
        behavior = parseDatafile_old(file);  % (entire workspace in this scope gets saved)
    end
    % save out results
    cd(getUserDataPath('pathType','Behavioral','session', Session,'run', Run));
    disp(['Saving behavioral data to ' pwd])
    save(SessionInfo.behaviorFileName(row),'behavior');
    cd(startPath)
end %End missingRunRows loop
end %End main function

% This version gets used for data logs older than January 20, 2020.
function t = parseDatafile_old(file)
disp('Now Parsing Data. Please wait...')

%find the beginning of trials
id_trial      = strfind(file,'trial');
rv_empty      = ~cellfun(@isempty,id_trial); %#ok<*STRCLFH>
id_trial_only = find(rv_empty>0);

%Initialization
nTrials = length(id_trial_only)-1;
cnt_attempt = 0;
wb = waitbar(0,'Parsing Behavior...');

for myTrial = 1:nTrials
    waitbar(myTrial/(length(id_trial_only)-1),wb,'Parsing Behavior...');
    
    % data from file, but only for the current trial
    file_current = file(id_trial_only(myTrial):id_trial_only(myTrial+1)-1);
    
    % only checking the first 30 data entries in this
    check_4_target = NaN(1,30);
    for myDummy           = 1:30
        Target                         = strfind(file_current(myDummy,:),'target');
        check_4_target(myDummy)        = ~cellfun(@isempty,Target);
    end
    
    % get the target position, e.g. [0.6000, 0.0]
    id_target_appear      =  find(check_4_target ==1);
    ReadTargetsPosition1  =  char(file_current(id_target_appear(1),:));
    left_parenth          =  strfind(ReadTargetsPosition1(1,:),'(');
    right_parenth         =  strfind(ReadTargetsPosition1(1,:),')');
    target_1              =  str2num(ReadTargetsPosition1(1,left_parenth+1:right_parenth-1)); %#ok<*ST2NM>
    
    % get the target prime position, e.g. [0.6000, 0.0]
    if length(id_target_appear)==2
        ReadTargetsPosition2  =  char(file_current(id_target_appear(2),:));
        left_parenth          =  strfind(ReadTargetsPosition2(1,:),'(');
        right_parenth         =  strfind(ReadTargetsPosition2(1,:),')');
        target_2              =  str2num(ReadTargetsPosition2(1,left_parenth+1:right_parenth-1));
    else
        target_2 = target_1;
    end
    
    % this may have been an attempted, but aborted, trial
    id_states           = strfind(file_current,'state');
    rm_empty_st         = ~cellfun(@isempty,id_states);
    id_states_only      = find(rm_empty_st>0);
    
    % set the state names
    states = {'initialfixation', 'fixationhold', 'cue', 'memory', ...
        'target_acquire', 'target_hold', 'reward'};
    
    if length(id_states_only) >= 5 % attempted (monkey got to memory)
        % set t(trial).trialID
        cnt_attempt                = cnt_attempt + 1;
        t(cnt_attempt).trialID     = cnt_attempt;
        
        % set t(trial) phase/state entry times
        t(cnt_attempt).trialstart        = TimeTrialStarts_old(file_current(1));   %#ok<*AGROW>
        for state = 1:length(id_states_only)-1
            eval(['t(cnt_attempt).' states{state} ' = TimeTrialStarts_old(file_current(id_states_only(state,:)));'])
        end
        % skip to the iti
        t(cnt_attempt).iti = TimeTrialStarts_old(file_current(id_states_only(length(id_states_only),:)));
        
        % set other attributes
        t(cnt_attempt).free_choice       = ~isequal(target_1,target_2);
        t(cnt_attempt).target            = sign(target_1(1));
        t(cnt_attempt).targetPos         = target_1;
        t(cnt_attempt).success           = length(id_states_only) > 7;
        t(cnt_attempt).lastState         = states{length(id_states_only)-1};
    end
end

% re order field for easier legibility
t = moveField(t, 'target_acquire', 'memory');
t = moveField(t, 'target_hold', 'target_acquire');
t = moveField(t, 'reward', 'target_hold');

close(wb);
end

%% sub function for parseDatafile_old. used for data logs older than January 20, 2020
function Start_Time = TimeTrialStarts_old(file_current)

Beh_data_current  = file_current;
id_time_start     = strfind(cell2mat(Beh_data_current(1,:)),' ');
tmp_t_s           = cell2mat(Beh_data_current(1,:));
Start_Time        = tmp_t_s(id_time_start(1)+1:id_time_start(2)-1);

end

%% begin parsing data. Used for data logs newer than january 20, 2020
function [t, EyePosMat, AnalogPosMat] = parseDatafile(file)

disp('Now Parsing Data. Please wait...')

% Pull out all eye and joystick position information
eyePositionInd = findCellArrayWithString(file(:,2),'eye:');
analogPositionInd = findCellArrayWithString(file(:,2),'analog:');
combinedPositionInd = eyePositionInd | analogPositionInd;
eyePositionInfo = file(eyePositionInd,:);
analogPositionInfo = file(analogPositionInd,:);


% Only run this portion of code if the user wants eye positions or analog positions though.
% It is time-intensive to run.
if nargout > 1
    %Parses time into matlab datenum format
    disp('Parsing eye positions');
    strDate = char(file(:,1));
    strDate = [strDate repmat(' ', size(strDate,1),1)];
    DateCell = textscan(strDate.', '%f-%f-%f %f:%f:%f,%f');
    DateMat = cell2mat(DateCell);
    DateVector = [DateMat(:,1:5) DateMat(:,6)+DateMat(:,7)/1000];
    parsedTime = datenum(DateVector);
    
    % Assemble record of eye positions [time xPos yPos]
    EyePosMat = zeros(length(eyePositionInfo),2);
    parfor i = 1:length(eyePositionInfo(:,2))
        EyePosMat(i,:) = parsePositionInput(eyePositionInfo{i,2});
    end
    EyePosMat = [parsedTime(eyePositionInd) EyePosMat];
    
    % Assemble record of joystick positions (first column is time, each
    % subsequent two rows are x,y for analog input
    
    % Currently hard-coding the 6 analog columns.
    numberOfAnalogInputs = 3;
    disp('Parsing analog positions');
    AnalogPosMat = zeros(length(analogPositionInfo), numberOfAnalogInputs * 2);
    parfor i = 1:length(analogPositionInfo(:,2))
        AnalogPosMat(i,:) = parsePositionInput(analogPositionInfo{i,2});
    end
    AnalogPosMat = [parsedTime(analogPositionInd) AnalogPosMat];
end

% Pull out trial info
trialInfo = file(~combinedPositionInd,:);
parfor i = 1:size(trialInfo,1)
    [trialInfo_messageType{i,1}, trialInfo_messageInfo{i,1}] = parseMessage(trialInfo{i,2});
end
trialInfo_parsed = horzcat(trialInfo(:,1), trialInfo_messageType, trialInfo_messageInfo);

% In this table:
% session: refers to number of distinct sessions within file itself
% trailNum: refers to trial number within the current session
% state: refers to state within current trial
% target: refers to target within current state
% parameter: refers to parameter name within current target
% value: refers to value of the parameter.

parameterTable = table('Size',[nnz(strcmp(trialInfo_parsed(:,2),'targetname')) 6],...
    'VariableTypes',{'double','double','string','string','string','string'},...
    'VariableNames',{'sessionNum','trialNum','state','target','parameter','value'});


session_startLine = find(strcmp(trialInfo_parsed(:,2),'task'));
parameter_counter = 1; % Counter for parameters in entire file
cnt_attempt = 0; %Counter for total number of successful trials in entire file (across multiple sessions within same file)
for session = 1:length(session_startLine)
    
    % data from file, but only for the current session
    if session ~= length(session_startLine)
        CurrentSession = trialInfo_parsed(session_startLine(session): session_startLine(session+1)-1,:);
    else
        CurrentSession = trialInfo_parsed(session_startLine(session):end,:);
    end
    
    % Fetch subject name
    subjectname = CurrentSession{strcmp(CurrentSession(:,2),'subject'),3};
    
    % Fetch task name
    taskname = CurrentSession{strcmp(CurrentSession(:,2),'task'),3};
    
    % Choose which task parse to use.
    % For new task types, need to create a new parser.
    switch taskname
        case 'fractal_value_overlap'
            overlap = true;
            [t,  parameterTable, parameter_counter, cnt_attempt] = ...
                MemoryMovementTask(CurrentSession, ...
                                  parameterTable, ...
                                  parameter_counter, ...
                                  cnt_attempt, ...
                                  session,...
                                  overlap);
        case 'fractal_value_memory'
            [t,  parameterTable, parameter_counter, cnt_attempt] = ...
                MemoryMovementTask(CurrentSession, ...
                                  parameterTable, ...
                                  parameter_counter, ...
                                  cnt_attempt, ...
                                  session);      
        case 'freechoice_memory_saccade'
            [t, parameterTable, parameter_counter, cnt_attempt] = ...
                MemoryMovementTask(CurrentSession, ...
                                parameterTable, ...
                                parameter_counter, ...
                                cnt_attempt, ...
                                session);    
        case 'instructed_mem_saccade_CorrectVersion'
            [t, parameterTable, parameter_counter, cnt_attempt] = ...
                MemoryMovementTask(CurrentSession, ...
                                parameterTable, ...
                                parameter_counter, ...
                                cnt_attempt, ...
                                session);
          case 'instructed_memory_movements'
            [t, parameterTable, parameter_counter, cnt_attempt] = ...
                MemoryMovementTask(CurrentSession, ...
                                parameterTable, ...
                                parameter_counter, ...
                                cnt_attempt, ...
                                session);      
        otherwise
            warning('This task type is not recognized and cannot be parsed.');      
    end %End switch statement for task parser
    
end % End session loop
end


%% subfunction to get eyeposition
function [position, units] = parsePositionInput(stringIn)
positionStartInd = strfind(stringIn,'(') + 1;
positionEndInd = strfind(stringIn, ')') - 1;
position = textscan(stringIn(positionStartInd:positionEndInd),'%f', 'Delimiter',',');
position = horzcat(position{:}');

unitStartInd = strfind(stringIn, '[') + 1;
unitEndInd = strfind(stringIn, ']') - 1;
%If less than 1, this format was broken and units cannot be parsed.
if unitStartInd > 1
    units = stringIn(unitStartInd:unitEndInd);
else
    units = [];
end
end

%% subfunction to get cell_ind for text_strings
function cell_ind = findCellArrayWithString(cellArrayIn,string)
temp_ind = strfind(cellArrayIn,string);
cell_ind = ~cellfun(@isempty,temp_ind);
end


%% subfunction to parse trial message
function [message_type, messageInfo] = parseMessage(stringIn)
colon_ind = strfind(stringIn,':');
message_type = stringIn(1:colon_ind(1)-1);

switch message_type
    case 'task'
        messageInfo = stringIn(colon_ind(1)+2:end-1);
    case 'subject'
        messageInfo = stringIn(colon_ind(1)+2:end-1);
    case 'trial'
        messageInfo = stringIn(colon_ind(1)+1:end-1);
    case 'choice'
        messageInfo = stringIn(colon_ind(1)+1:end-1);
    case 'RWD magnitude'
        messageInfo = stringIn(colon_ind(1)+1:end-1);
    case 'blocknum'
        messageInfo = stringIn(colon_ind(1)+1:end-1);
    case 'effector'
        % Outputs array as {effectorType windowConstraint}
        window_constraint_ind = strfind(stringIn,'window_constraint:');
        messageInfo.effector = stringIn(colon_ind(1)+1:window_constraint_ind-2);
        messageInfo.windowconstraint = stringIn(window_constraint_ind+18:end-1);
    case 'imageID'
        messageInfo = stringIn(colon_ind(1)+2:end-1); %Image ID
    case 'target'
        % outputs array as {targetname {xPos yPos} unit}
        target_type = stringIn(colon_ind(1)+1:colon_ind(2)-1);
        [targetPos, units] = parsePositionInput(stringIn);
        messageInfo.targetName = target_type;
        messageInfo.targetPos = targetPos;
        messageInfo.targetUnits = units;
    case 'state'
        messageInfo = stringIn(colon_ind(1)+1:end-1); %state
    case 'targetname'
        % Outputs target parameters as {targetname parametername value}
        parameter_ind = strfind(stringIn,'parameter:');
        value_ind = strfind(stringIn,'value:');
        targetname = stringIn(colon_ind(1)+2:parameter_ind-2);
        parametername = stringIn(parameter_ind + 11:value_ind-2);
        value = stringIn(value_ind+7:end-1);
        messageInfo.targetname = targetname;
        messageInfo.parameter = parametername;
        messageInfo.value = value;
    
    otherwise
        messageInfo = [];
        
end

end

