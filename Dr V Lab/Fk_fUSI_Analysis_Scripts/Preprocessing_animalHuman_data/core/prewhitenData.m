% Script that prewhitens data for each session/run
% 
% downsamplefactor = how much to downsample each image.
%
% Loosely based upon code from Vasileios


% User settings
downsamplefactor=40;

% load Session Info
load('SessionInfo.mat','SessionInfo');


%% get the sessions & runs desired

% create the list of available sessions/runs
runStrings = cell(size(SessionInfo,1),1);
counter=0;
for i = 1:size(SessionInfo,1) 
    Session = SessionInfo.Session(i);
    Run = SessionInfo.Run(i);
    runStrings{i} = ['Session ' num2str(Session) ', Run ' num2str(Run)];
end

% use a GUI to select which ones you want
tf = 0;
while ~tf
    [indx,tf] = listdlg('PromptString','Select Session/Run to load:',...
        'SelectionMode','multiple',...
        'ListString',runStrings);
end
sessionRunList = [SessionInfo.Session(indx) SessionInfo.Run(indx)];

%Don't use these sessions/runs
dontuse = [3 1; 3 2; 3 3; 4 1; 4 2; 5 1; 5 2; 19 1;1 3];
questionuse = [14 1; 12 1; 12 2;];


load('BadRunList.mat','badRuns');
dontuseInd = ismember(sessionRunList,badRuns,'rows');
sessionRunList = sessionRunList(~dontuseInd,:);

%% Load ROIs
AnatomicalPolygons_fileLocation = getUserDataPath('pathType','AnatomicalPolygons');
AnatomicalPolygons_vertices = loadAnatomicalPolygons(AnatomicalPolygons_fileLocation,sessionRunList);


%% Load Doppler Data and pre-whiten every voxel
for i=1:size(sessionRunList,1)
    session = sessionRunList(i,1);
    run = sessionRunList(i,2);
    currentRow = intersect(find(SessionInfo.Session==session),...
        find(SessionInfo.Run==run));
    [iDop, LInds, RInds, behavior, coreParams, wallpaper] = concatenateRuns([session run],true);
    iDopP = preProcess(iDop, 'timeGain', true, 'diskFilter', true, 'zScore', false,'downsample',1/downsamplefactor);
    [yPixels, xPixels, nWindows, nTrials] = size(iDopP);
    
    %Pre-whiten data
    iDopRes = reshape(iDopP,size(iDopP,1),size(iDopP,2),[]);
    
    p = gcp('nocreate');
    if isempty(p)
        parpool;
    end
    ppm = ParforProgMon('Parfor Loop Progress', size(iDopRes,1));
    
    tic
    initial_arima_model = arima(10,1,1);
    iDopInov = NaN(size(iDopRes,1),size(iDopRes,2),size(iDopRes,3));
    n = size(iDopRes,2);
    parfor ii = 1:size(iDopRes,1)
        iDopInov_temp=[];
        for jj = 1:n
            arima_model = estimate(initial_arima_model,squeeze(iDopRes(ii,jj,:)),'Display','off');
            [E,~] = infer(arima_model,squeeze(iDopRes(ii,jj,:)));
            iDopInov_temp(jj,:) = E;
        end
        iDopInov(ii,:,:)=iDopInov_temp;
        ppm.increment();
    end
    toc
    iDopInov = reshape(iDopInov,yPixels,xPixels,nWindows,nTrials);
    foldername=sprintf('%s%dX downsample',getUserDataPath('pathType','PrewhitenedData'),downsamplefactor);
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end
    save(fullfile(foldername,sprintf('Session%d_Run%d',foldername,session, run)),'iDopInov','downsamplefactor');
end
