% baseline_drift loads data using 'loadSingleRunData and then uses a basic
% linear fit to characterize each voxel's drift as a function of successful
% trial (i.e. through time). The results are then plotted as a percent
% drift from the beginning of the trial to the end (subplot 1), r^2 of the
% fit (confidence of each voxel's fit), and finally a reference image of
% mean fUS activation (sanity check). 
% 
% This script can be run without passing any data manually. 

%% load doppler data 
% loadDopplerData(whos)
% [yPixels, xPixels, nWindows, nTrials] = size(iDop);
% nDataPts = nWindows*nTrials;

iDop = permute(fL5_preMC, [2 1 3 4]);
[yPixels, xPixels, nWindows, nTrials] = size(iDop);
nDataPts = nWindows*nTrials;

%% get the r squared and percent change maps, voxel by voxel (linear fits)
[rsq, percChange] = deal(zeros(yPixels,xPixels));

thePool = gcp;
H = ParforProgMon('Measuring drift. Please wait: ', yPixels, 1);
parfor y = 1:yPixels
    for x = 1:xPixels
   
        voxelData = squeeze(iDop(y,x,:,:));         
        voxelData = reshape(voxelData,[],1);
        coeffs = polyfit((1:nDataPts)', voxelData, 1);        
        yFit = polyval(coeffs, 1:nDataPts);
        percChange(y,x) = (yFit(end)-yFit(1))/yFit(1)*100;
        
        % Plot the fitted line (debugging only)
%         clf; hold on;
%         plot(voxelData)
%         plot(1:nDataPts, yFit, 'r-', 'LineWidth', 3);
%         xlabel('successful trial'); ylabel('baseline power'); pause(0.01);
        
        yresid = voxelData'-yFit;    % compute the residual values
        SSresid = sum(yresid.^2);   % residual sum of squares 
        % Total sum of squares = variance multiplied by number of observations minus 1:
        SStotal = (nDataPts-1) * var(voxelData);
        rsq(y,x) = 1 - SSresid/SStotal;  % compute R2 
        
    end
    H.increment();
end

%% plot the results
X_img_mm = (0:size(iDop,2)-1)*0.1/4;
Z_img_mm = (0:size(iDop,1)-1)*0.1/4;

% plot the percent change over 62 trials
figure('position',[50 400 1850 500])
subplot(131)
imagesc(X_img_mm, Z_img_mm, percChange, [-300 300]); colorbar
xlabel('mm'), ylabel('mm'), title(['Percent Change over ' num2str(nDataPts) ' trials'])
daspect([max(X_img_mm) max(Z_img_mm) 1])

% plot the fit r^2 for each voxel
subplot(132)
imagesc(X_img_mm, Z_img_mm, rsq, [0 1]); colorbar
xlabel('mm'), ylabel('mm'), title('Line fit r^2')
daspect([max(X_img_mm) max(Z_img_mm) 1])

% plot the basic vascular map (for reference)
subplot(133)
Ibg = squeeze(mean(iDop(:,:,1,:),4)); % background vascular US image
Ibg = (Ibg - min(Ibg(:)))./(max(Ibg(:))-min(Ibg(:))); % normalize to [0 1]
Ibg = nthroot(Ibg,5);                    % nonlinear scale between [0 1]
imagesc(X_img_mm, Z_img_mm, Ibg); 
xlabel('mm'), ylabel('mm'), title('Average vascular map over session')
daspect([max(X_img_mm) max(Z_img_mm) 1])

% put some info down
if coreParams.motionCorrection
    suptitle(['S' num2str(Session) 'R' num2str(Run) ', ''' coreParams.method ''' correction'])
else
    suptitle(['S' num2str(Session) 'R' num2str(Run) ', no motion correction'])
end

% plot the percent change mulitplied by r^2
% subplot(133)
% imagesc(X_img_mm, Z_img_mm, percChange.*rsq, [-300 300]); colorbar
% xlabel('mm'), ylabel('mm'), title('(percent change) \times r^2')
% daspect([max(X_img_mm) max(Z_img_mm) 1])

colormap hot