%% load some sample doppler data
loadDopplerData(whos)
[yPixels, xPixels, nWindows, nTrials] = size(iDop);

%% display a video across time
response = input('Do you want to see a video of the doppler data? (y/n)  ','s');

% continuousTimeDop =  reshape(iDop_lc,[yPixels, xPixels, nWindows*nTrials]);
if strcmp(response,'y')
    figure()
    for trial = 1:nTrials
        for window = 1:nWindows
            % get the current frame
            frame = nonLinearScaling(iDop(:,:,window,trial));        

            % plot it
            imagesc(frame), colormap hot
            xticklabels({''}), yticklabels({''})
            title(['Trial = ' num2str(trial) '/' num2str(nTrials)])
            pause(0.01);
        end
    end
end

%% try a basic motion correction on fUS data using both nonrigid & imregister

first = iDop(:,:,1,1);
fixed = iDop(:,:,ceil(nWindows/2),ceil(nTrials/2));
moving = iDop(:,:,nWindows,nTrials);

fprintf('Registering image(s)...')
% Perform the registration using imregister
[optimizer, metric] = imregconfig('monomodal');
% optimizer.InitialRadius = 0.009;      % only for use in multimodal
% optimizer.Epsilon = 1.5e-4;           % only for use in multimodal
% optimizer.GrowthFactor = 1.01;        % only for use in multimodal
optimizer.MaximumIterations = 300;
movingRegistered = imregister(moving, fixed, 'affine', optimizer, metric);

% Perform the registration using rigid transform
Options.Similarity = 'sd';
Options.Registration = 'rigid';
[~, Grid, Spacing] = image_registration(moving, fixed, Options);
movingRegisteredRigid = bspline_transform(Grid, moving, Spacing);

% Perform the registration using normcorre
iDopNormcorre = normcorre_doppler(iDop);
movingRegisteredNormcorre = iDopNormcorre(:,:,nWindows,nTrials);

%% plot the images
close all, figure('position',[100 50 1600 600])
subplot(221), imshowpair(nonLinearScaling(fixed), nonLinearScaling(moving),'Scaling','joint'), title('un-registered overlay')
subplot(222), imshowpair(nonLinearScaling(fixed), nonLinearScaling(movingRegisteredRigid),'Scaling','joint'), title('registered by rigid')
subplot(223), imshowpair(nonLinearScaling(fixed), nonLinearScaling(movingRegistered),'Scaling','joint'), title('registered by imregister')
subplot(224), imshowpair(nonLinearScaling(fixed), nonLinearScaling(movingRegisteredNormcorre),'Scaling','joint'), title('registered by normcorre')

%% function to get nonlinear scaling for easier viewing of fUS images
function imOut = nonLinearScaling(imIn)
        imIn = (imIn - min(imIn(:)))./(max(imIn(:))-min(imIn(:))); % normalise background data [0 1]
        imOut = nthroot(imIn,4);           % nonlinear scale between [0 1]       
end
