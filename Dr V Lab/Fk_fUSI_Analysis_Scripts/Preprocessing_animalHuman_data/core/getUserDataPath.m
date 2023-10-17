function path = getUserDataPath(varargin)
% path = getUserDataPath(varargin)
%
% getUserDataPath automatically determines directory to that of the data
% type specified, e.g. 'Behavioral' or 'fUS', for the Session and Run of
% interest.
%
% You can also invoke pathType 'fast' to load from a secondary
% (probably local) directory that is faster, e.g. a local SSD. Note that
% these 'fast' directories will only include the S##R##.mat doppler data
% and will not have access to raw data!
%
% %%%%%%%%%%%%%%%%%%%  VARARGIN:
% pathType accepts any of following inputs
%
% 'default': root path, e.g. Z:\fUS\
% 'Processed': processed data path, e.g. Z:\fUS\Processed Data\
% 'Raw': raw data path, e.g. Z:\fUS\Raw Data\Raw Data\
% 'deepfus': deep fUS, e.g. Z:\fUS\deep fUS\
% 'models': trained ML models, e.g. Z:\fUS\deep fUS\pretrained_models\
%
% these options append the processed data path:
% 'ExtractedData', 'PrewhitenedData', 'AnatomicalPolygons', 'SulcusMap'
% 'Figure' or 'Output' both go to  path to output folder, e.g. Z:\fUS\Processed Data\output
% 'ROIs': path to ROIs
% 'anatomicalpolygons'
% 'sulcusmap'
% 'prewhiteneddata'
%
% paths that refer to SessionInfo as a lookup table for session/run specific paths:
% 'Behavioral': path to raw behavioral data, e.g. Z:\fUS\Raw Data\20201231\
% 'fUS': path to raw (.bin) fUS data, e.g. Z:\fUS\Raw Data\20201231\Acq_121535
% * note that you must pass session & run info with these options!
%
% session: session number (positive integer)
%
% run: run number (positive integer)
%
% speed: 'default' or 'fast' (string)
%
% model: 'nhp' or 'human' (string)
%
% %%%%%%%%%%%%%%%%%%%  VARARGOUT:
% path: path to the desired data
%
% %%%%%%%%%%%%%%%%%%%  EXAMPLE:
% This function uses input parser. Example usage:
%
% path =
% getUserDataPath('pathType','Behavior','speed','fast','session',2,'run',1)
% would go to the path with the behavioral data for session 2 run 1

%% Variable argument input parser
p = inputParser;
p.addParameter('pathType','default', @ischar)
p.addParameter('session','', @isnumeric)
p.addParameter('run','', @isnumeric)
p.addParameter('speed','default',@ischar)
p.addParameter('model','nhp',@ischar)
p.parse(varargin{:});
VariableInputs = p.Results;

% we don't care about case-sensitivity
VariableInputs.pathType = lower(VariableInputs.pathType); 
VariableInputs.speed = lower(VariableInputs.speed);
VariableInputs.model = lower(VariableInputs.model);

%% load session data paths
try
    if strcmp(VariableInputs.model,'human')
        load('SessionInfoTBI.mat','SessionInfo')
    else
        load('SessionInfo.mat','SessionInfo')
    end
catch
    if nargin == 0
        warning('SessionInfo.mat (library of data paths) not found!')
    elseif nargin >= 1
        error('SessionInfo.mat (library of data paths) not found!')
    end
end

%% create end of path
switch VariableInputs.pathType
    case 'processed'
        if strcmp(VariableInputs.model,'human')
            appendPath = 'Processed TBI Data';
        else
            appendPath = 'Processed Data';
        end
    case 'raw'
        appendPath = 'Raw Data';
    case {'behavioral','fus'}
        % make sure session/run is defined & get row index from SessionInfo
        session = VariableInputs.session;
        run = VariableInputs.run;
        if isempty(session) || isempty(run)
            error('Did not specify session or run to load.');
        else
            rowOfInterest = intersect(  find(SessionInfo.Session==session),...
                find(SessionInfo.Run==run));
        end
        
        % check for multiple matching entries
        if length(rowOfInterest) > 1
            warning(['More than one Session/Run entry exists!' ...
                'Are you sure you haven''t processed this data before?'])
            rowOfInterest = rowOfInterest(length(rowOfInterest));
        end
        
        % check for no matching enries or invalid row numbers
        if isempty(rowOfInterest) || ...
                rowOfInterest > size(SessionInfo,1) || ...
                rowOfInterest <= 0
            error('Session/Run combination is invalid!')
        end
        
        % creating path of type: 'Raw Data\Behavior Path\fUS path' for NHP
        % creating path of type: 'Raw Data\fUS path' for human
        appendPath = 'Raw Data';
        if strcmp(VariableInputs.pathType,'behavioral')
            appendPath2 = SessionInfo.Behavioral{rowOfInterest};
            appendPath = fullfile(appendPath,appendPath2);
        elseif strcmp(VariableInputs.pathType,'fus') && strcmp(VariableInputs.model,'nhp')
            appendPath2 = SessionInfo.Behavioral{rowOfInterest};
            appendPath3 = SessionInfo.fUS{rowOfInterest};
            appendPath = fullfile(appendPath,appendPath2,appendPath3);
        elseif strcmp(VariableInputs.pathType,'fus') && strcmp(VariableInputs.model,'human')            
            appendPath2 = SessionInfo.fUS{rowOfInterest};
            appendPath = fullfile(appendPath,appendPath2);
        end    
    case 'rois'
        appendPath = fullfile('Processed Data','ROIs');
    case 'anatomicalpolygons'
        appendPath = fullfile('Processed Data','Anatomical_Polygons');
    case 'sulcusmap'
        appendPath = fullfile('Processed Data','SulcusMaps');
    case 'prewhiteneddata'
        appendPath = fullfile('Processed Data','Prewhitened data');
    case {'figure','output'}
        if strcmp(VariableInputs.model,'human')
            appendPath = fullfile('Processed TBI Data','output');
        else
            appendPath = fullfile('Processed Data','output');
        end
    case 'extracteddata'
        appendPath = fullfile('Processed Data','Extracted Data');
    case 'default'
        appendPath = '';
    case 'models'
        appendPath = 'Deep fUS\pretrained_models';
    case 'deepfus'
        appendPath = 'Deep fUS\';
    otherwise
        error('path type problem: path not set!');
end

%% get user & computer names
user = getenv('username');
computer = getenv('computername');

% some older Unix OS versions will return empty arrays for getenv...
if isempty(user) || isempty(computer)
    try
        user = char(java.lang.System.getProperty('user.name'));
        %Not currently used, but leaving in case this is useful later.
        computer = char(java.net.InetAddress.getLocalHost.getHostName);
    catch
        error('Could not get username')
    end
end

%% set path according to user & computer

% sumner @ caltech
if strcmp(user,'Sumner') && ~strcmp(VariableInputs.speed,'fast')
    if strcmp(VariableInputs.speed,'fast')
        basePath = 'F:\';
    else
        basePath = 'Z:\fUS';
    end
    
% sumner @ Stardust    
elseif strcmp(user,'sumner') && contains(computer,'Stardust')
        if  strcmp (VariableInputs.speed,'fast')
            basePath = '/Users/sumner/Documents/data';
        else
            basePath = '/Volumes/fUS_2020';            
        end
        
% whitney     
elseif strcmp(user,'wsgriggs') 
    if strcmp(VariableInputs.speed,'fast')
        basePath = 'C:\Users\wsgriggs\Box\Andersen Lab\Data\fUS';
    else
        basePath = 'Z:\fUS';
    end
elseif strcmp(user,'WhitneyGriggs')
    basePath = '/Volumes/raid01/fUS';
    
% claire    
elseif strcmp(user,'Claire')
    basePath = 'D:\fUS human\Data\';
else
    error('path has not been set for this user!')
end


%% create full path 

path = fullfile(basePath,appendPath);

% make sure our file separators are correct
path = strrep(path,'/',filesep);
path = strrep(path,'\',filesep);

end