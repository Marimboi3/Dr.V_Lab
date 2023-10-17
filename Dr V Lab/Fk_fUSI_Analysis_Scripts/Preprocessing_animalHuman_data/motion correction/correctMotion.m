function [dopOut, coreParams,template] = correctMotion(dopIn, coreParams, verbose, template)
% this function uses rigid body transforms to correct errant motion in the
% fUS dopppler sequence. It does this on a single-trial basis, i.e. the
% data passed to this function is of the size X x Y x nWindows x nTrials.
% nWindows is the number of fUS activation images in a single trial, e.g.
% ~20 images. 
% 
% inputs:
%
% dopIn: xPixels x yPixels x nWindowsInEpoch x nTrials
% coreParams: the coreparams structure used in createDoppler (and saved in
%             S##R##.mat doppler data files). note: coreParams should
%             contain coreParams.method of a valid choice: 'imregister',
%             'rigid', or 'normcorre'. Otherwise a default is chosen.
%
% output:
%
% dopOut: array xPixels x yPixels x nWindowsInEpoch x nTrials  (corrected)
% coreParams: will now contain the reference frame and trial 

%% check the registration method
if ~isfield(coreParams,'method')
    warning('defaulting to normcorre method!')
    method = 'normcorre';
else 
    method = coreParams.method;
end

if ~strcmp(method,'imregister') && ~strcmp(method,'rigid') && ~strcmp(method,'normcorre')    
end

if ~exist('verbose','var')
    verbose = false;
end

if ~exist('template','var')
    useTemplate = false;
else
    useTemplate = true;
end

%% select reference frame and window
[~, ~, nWindows, nTrials] = size(dopIn);
% refTrial: the trial number to pull the reference (stationary) frame from
% refFrame: the frame/window within the trial  
if strcmp(method,'imregister') || strcmp(method,'rigid')
    refTrial = ceil(nTrials/4);
    refFrame = 1;
    fixedFrame = dopIn(:,:,refFrame,refTrial);
elseif strcmp(method,'normcorre')
    refTrial = 0; 
    refFrame = 0;
else
    error('invalid method chosen.')
end

fprintf('Reference: trial %i/%i, frame %i/%i\n', refTrial, size(dopIn,4), refFrame, size(dopIn,3))
% save out updated coreParams structure
coreParams.motionCorrection_refTrial = refTrial;
coreParams.motionCorrection_refFrame = refFrame;

%% choose the reference frame & initialize data
dopOut = zeros(size(dopIn));

%% imregister method
if strcmp(method,'imregister')
    [optimizer, metric] = imregconfig('monomodal');
    
%     optimizer.InitialRadius = 0.009; % only for use in multimodal
%     optimizer.Epsilon = 1.5e-4; % only for use in multimodal
%     optimizer.GrowthFactor = 1.01; % only for use in multimodal
    optimizer.MaximumIterations = 300;
    
    % for each frame in the sequence
    wb = waitbar(0,'Correction motion using imregister...');
    for trial = 1:nTrials
        for window = 1:nWindows
            
            % select the current frame to work on
            movingFrame = dopIn(:,:,nWindows,nTrials);
            
            % register the image
            dopOut(:,:,window,trial) = imregister(movingFrame, fixedFrame, 'affine', optimizer, metric);
            
            % update the waitbar
            waitbar(trial/size(dopIn,4));
        end
    end
    close(wb)
end

%% rigid method
if strcmp(method,'rigid')
    Options.Similarity = 'sd';
    Options.Registration = 'rigid'; %#ok<*STRNU>
    
    % for each frame in the sequence...
    wb = waitbar(0,'Correction motion using rigid transform...');
    for trial = 1:nTrials
        for window = 1:nWindows
            
            % select the current frame to work on
            movingFrame = dopIn(:,:,window,trial);
            
            % register the image (wrapped evalc to prevent massive cmd prints)
            evalc('[dopOut(:,:,window,trial), Grid, Spacing] = image_registration(movingFrame, fixedFrame, Options);');
            waitbar(trial/size(dopIn,4));
            
            % transform the result & save to dopOut
            dopOut(:,:,window,trial) = bspline_transform(Grid, movingFrame, Spacing);
        end
    end
    close(wb)
end

%% normcorre method
if strcmp(method,'normcorre')
    disp('Correcting motion using normcorre')
    if useTemplate
        [dopOut,template] = normcorre_doppler(dopIn, verbose,template);
    else
        [dopOut,template] = normcorre_doppler(dopIn, verbose);
    end
end

end % function

