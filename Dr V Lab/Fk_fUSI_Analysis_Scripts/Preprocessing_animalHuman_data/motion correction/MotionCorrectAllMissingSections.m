

clear; close all; format compact
startDir = pwd;


%% parameters & useful vars
coreParams.motionCorrection = true;
coreParams.missedTrials = false;
% set the type of motion correction used (rigid, imregister, normcorre)
coreParams.method = 'normcorre';   
coreParams.genDate = date;

%% run through the entire SessionInfo file to find which data has not been motion corrected
load SessionInfo.mat
nEntries = size(SessionInfo,1);
fprintf('SessionInfo.mat loaded successfully. %i entries found.\n\n', nEntries)    

missingRunRows = [];
% move to the processed data folder
cd(getUserDataPath('pathType','Processed')); 

% Check to see which motion corrected files are missing
for row = 1:size(SessionInfo,1)        
    Session = SessionInfo.Session(row);
    Run = SessionInfo.Run(row);
    
    dopplerFileName = makeFileString(Session, Run, coreParams);

    if ~exist(dopplerFileName,'file') 
        missingRunRows = [missingRunRows row]; %#ok<AGROW>
        fprintf('Session %i Run %i %s is missing.\n\n',...
            Session, Run, dopplerFileName)
    else
        fprintf('Session %i Run %i %s found!\n',...
            Session, Run, dopplerFileName)
    end    
end

fprintf('%i out of %i doppler files found!\n',...
    nEntries - length(missingRunRows), nEntries)
if ~isempty(missingRunRows)
    uiwait(msgbox('Click OK to generate missing files (ctrl+c to quit)'),120);
end

%% Load in non-motion corrected Doppler data, motion correct, and save data
for row = missingRunRows
    Session = SessionInfo.Session(row);
    Run = SessionInfo.Run(row);
    
    fprintf('Motion correcting Session %d Run %d \n',Session, Run);
    
    LoadParams.motionCorrection = false; % For loading non-motion corrected data
    LoadParams.missedTrials = false;
    dopplerFileName = makeFileString(Session, Run, LoadParams);
    load(fullfile(getUserDataPath('pathType','Processed'), dopplerFileName))

    % Setting necessary core params to desired values for motion correction
    coreParams.motionCorrection = true;
    % set the type of motion correction used (rigid, imregister, normcorre)
    coreParams.method = 'normcorre';
    coreParams.genDate = date;
    
    [iDop, coreParams] = correctMotion(iDop, coreParams); 
    
    fileString = makeFileString(Session, Run, coreParams);
    cd(getUserDataPath('pathType','Processed'));
    fprintf('Saving %s to %s\n', fileString, pwd)
    save(fileString, 'iDop', 'coreParams', 'UF', 'behavior','LInds','RInds', 'wallpaper', '-v7.3')    
    
end

function fileString = makeFileString(Session, Run, coreParams)
% Same function as from createDopplerImages.m
    baseString = append('doppler_S', num2str(Session), '_R', num2str(Run));
    if coreParams.missedTrials
        baseString = append(baseString, '_missedTrials');
    end
    if coreParams.motionCorrection        
        baseString = append(baseString, '+', coreParams.method);
    end       
    fileString = append(baseString,'.mat');
end