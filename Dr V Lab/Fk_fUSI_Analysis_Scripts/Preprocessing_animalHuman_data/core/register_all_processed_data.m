% register_all_processed_data is a simple script to load any Session/Run
% combinations that do not already have a pre-registered partner in the
% database and create that partner.

function register_all_processed_data()
clear; close all

%% look for missing datasets
load('SessionInfo.mat','SessionInfo');
missingRunRows = [];
startPath = pwd;
cd(getUserDataPath('pathType','Processed','speed','fast'));
for row = 1:size(SessionInfo,1)
    Session = SessionInfo.Session(row);
    Run = SessionInfo.Run(row);
    
    % if we can't find the doppler file where it belongs...
    dopplerFileName = ['doppler_S' num2str(Session) '_R' num2str(Run) '+normcorre.mat'];
    if ~exist(dopplerFileName,'file')
        missingRunRows = [missingRunRows row]; %#ok<AGROW>
        fprintf('Session %i Run %i %s is missing.\n\n',...
            Session, Run, dopplerFileName)
    else
        fprintf('Session %i Run %i %s found!\n',...
            Session, Run, dopplerFileName)
    end
end

%% create missing doppler files

for row = missingRunRows
    % load the unregistered data
    Session = SessionInfo.Session(row);
    Run = SessionInfo.Run(row);
    fprintf('\nRegistering Session %i Run %i:\n', Session, Run);
    dopplerFileName = ['doppler_S' num2str(Session) '_R' num2str(Run) '.mat'];
    load(dopplerFileName); %#ok<LOAD>
    
    % create pre-registration wallpaper
    if ~exist('wallpaper','var')
        wallpaper = nthroot(squeeze(mean(mean(iDop,3),4)),3);
        lowerCutoff = quantile(wallpaper(:),0.01);
        wallpaper(wallpaper<lowerCutoff) = lowerCutoff;
    end
    
    % motion correct the data
    coreParams.method = 'normcorre';
    coreParams.motionCorrection = true;
    [iDop, coreParams] = correctMotion(iDop, coreParams);
    
    % save the data
    saveString = strrep(dopplerFileName, '.mat', '+normcorre.mat');
    save(saveString, 'iDop', 'coreParams', 'UF', 'behavior','LInds','RInds','wallpaper', '-v7.3') %#ok<USENS>
end
cd(startPath)

end