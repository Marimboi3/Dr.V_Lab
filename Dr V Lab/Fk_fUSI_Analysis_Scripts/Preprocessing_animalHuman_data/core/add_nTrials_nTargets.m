clear; close all
startPath = pwd; % note: must be ran from base path! 
load SessionInfo.mat
disp(SessionInfo)
startPath = cd(getUserDataPath('pathType','Processed'));
for row = 60:size(SessionInfo,1)
    session = SessionInfo.Session(row);
    run = SessionInfo.Run(row);
    disp(['Loading Session ' num2str(session) ', run ' num2str(run)])
    filename = strcat('doppler_S', num2str(session), '_R', num2str(run), '+normcorre.mat');
    load(filename);
    
    nTrials = size(iDop,4);    
    
    targetPos = vertcat(behavior.targetPos);
    [angle,distance] = cart2pol(targetPos(:,2),targetPos(:,1)); % convert to polar coordinates
    angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.
    nTargets = length(unique(angle));    
    
    SessionInfo.nTrials(row) = nTrials;
    SessionInfo.nTargets(row) = nTargets;
end
disp(SessionInfo)
uiwait(msgbox('Are you happy?'))
save('SessionInfo.mat','SessionInfo');
cd(startPath)
save('SessionInfo.mat','SessionInfo');
