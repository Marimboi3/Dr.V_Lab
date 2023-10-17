function dopOut = preProcess(dopIn, varargin)
% dopOut = preProcess(dopIn, vargin)
%
% preProcess will perform basic pre-processing steps on the raw doppler
% data of the form: [yPixels x xPixels x timeWindows x trials]
% Good practice would require that all doppler data is passed together,
% regardless of target. For example, pass iDop, but not dopIn.
%
% important: preProcess is usually called immediately after loading a
% session's data using loadSessionData
%
% optional vargin:
% 'interp',     bool 
% 'uninterp',   bool
% 'timeGain',   bool
% 'zscore',     bool
% 'diskfilter', bool
% 'downsample', a 1x2 int array, depth x width, e.g. [200 200]
% 'preWhiten',  bool
% 'detrend', int array of trial time baseline indices, e.g. 5:8
% 'verbose',    bool
%
% vargout
% dopOut is of the same form as dopIn, but with preprocessing steps

%% handling varargin
p = inputParser;
p.addOptional('interp',false,@islogical)
p.addOptional('uninterp',false,@islogical)
p.addOptional('timeGain',false,@islogical)
p.addOptional('zScore',false,@islogical)
p.addOptional('diskFilter',false,@islogical)
p.addOptional('downsample',false)
p.addOptional('preWhiten',false,@islogical)
p.addOptional('crop',false)
p.addOptional('detrend',false)
p.addOptional('verbose',false,@islogical)
p.parse(varargin{:});
sets = p.Results;
if sets.verbose
    disp(sets)
end


%% main function is here (function calls)
if sets.interp, dopIn = interDop(dopIn); end
if sets.uninterp, dopIn = uninterpDop(dopIn); end
if sets.timeGain, dopIn = timeGainCompensation(dopIn,sets); end
if sets.zScore, dopIn = dopZ(dopIn); end
if sets.diskFilter, dopIn = diskFilter(dopIn); end
if sets.downsample, dopIn = imresize(dopIn, sets.downsample); end
if sets.preWhiten, dopIn = white(dopIn); end
if sets.detrend, dopIn = detrendDrift(dopIn, sets.detrend); end
if sets.crop, dopIn = cropEdges(dopIn, sets.crop); end

dopOut = dopIn;

%% interpolation
    function dataOut = interDop(dataIn)
        [n_depth, n_width, n_windows, n_trials] = size(dataIn);
        dataOut = nan(4*n_depth-3, 4*n_width-3, n_windows, n_trials);
        for t = 1:n_trials
            for w = 1:n_windows
                dataOut(:,:,w,t) = interp2(dataIn(:,:,w,t),2);
            end
        end
    end

%% uninterpolation
    function dataOut = uninterpDop(dataIn)
        % note that this uninterpolation only works to undo the use of
        % interp2(data,2), as was usually done in the create Doppler core
        % prior to November 2020.
        [n_depth, n_width, n_windows, n_trials] = size(dataIn);        
        dataOut = nan((n_depth+3)/4, (n_width+3)/4, n_windows, n_trials);
        for t = 1:n_trials
            for w = 1:n_windows
                dataOut(:,:,w,t) = dataIn(1:4:end,1:4:end,w,t);
            end
        end
    end

%% disk filter
    function dataOut = diskFilter(dataIn)      
        [~, ~, n_windows, n_trials] = size(dataIn); 
        dataOut = NaN(size(dataIn));
        h = fspecial('disk',2); % disk size defined here
        
        % can't filter2 n-d arrays -> ugly nested for-loop (sorry)
        % note: NOT faster with parfor overhead (I tried)
        for window = 1:n_windows
            for trial = 1:n_trials
                dataOut(:, :, window, trial) = ...
                    filter2(h, squeeze(dataIn(:, :, window, trial)));
            end
        end        
    end

%% time gain compensation
    function dataOut = timeGainCompensation(dataIn, res)
        % getting useful dimensions        
        [n_depth, n_width, n_windows, n_trials] = size(dataIn); 
        depthInds = (1:n_depth)';    
        indicesForModeling = round(0.1*n_depth):n_depth; % discards supra-dural info
        
        % get mean signal (empirical)
        dataIn_flat = reshape(dataIn, n_depth, n_width, []);
        depthMeanSignal = squeeze(mean(mean(dataIn_flat,2),3)); % 1D depth array (actual)
        
        % create the exponential model for depth attenuation/time-gain-comp        
        f = fit(depthInds(indicesForModeling) , double(depthMeanSignal(indicesForModeling)), 'exp1');
        depthMean = f.a*exp(f.b*depthInds);
        
        % plot the result
        if res.verbose
            figure();
            plot(depthInds, depthMean, 'k', 'lineWidth', 2);
            hold on; plot(depthInds, depthMeanSignal, '.r')
            legend('depth attenuation model','actual mean activation')
            xlabel('depth (pixels)'); ylabel('raw signal [a.u.]')
        end
        
        % normalization happens here
        dataOut = dataIn./repmat(depthMean, 1, n_width, n_windows, n_trials);               
        dataOut = dataOut*100;  % converts to percentage of mean signal                      
    end

%% z score each voxel
    function dataOut = dopZ(dataIn)        
        [n_depth, n_width, n_windows, n_trials] = size(dataIn); 
                      
        % flatten the data into a voxel-column for zscoring
        dataOut = permute(dataIn, [2, 1, 4, 3]);% [n_width, n_depth, n_trials, n_windows]);
        dataOut = reshape(dataOut, [n_depth*n_width, n_windows*n_trials]);        
        % perform the zscore 
        dataOut = zscore(dataOut')';
        % reshape back to original data dims
        dataOut = reshape(dataOut, [n_width, n_depth, n_trials, n_windows]);        
        dataOut = permute(dataOut, [2, 1, 4, 3]);                
    end

%% pre-whiten the data
    function dataOut = white(dataIn)
        [n_depth, n_width, n_windows, n_trials] = size(dataIn); 
        
        % flatten the data across time
        dataOut = permute(dataIn, [2, 1, 4, 3]);% [n_width, n_depth, n_trials, n_windows]);
        dataOut = reshape(dataOut, [n_depth*n_width, n_windows*n_trials]);               
        
        % create an arima model 
        mdl = arima(10,1,1);      
        imageMean = transpose(double(squeeze(dataOut(ceil(end/2),:)))); % uses center-most px
        EstMdl = estimate(mdl, imageMean, 'Display', 'off');
        % correct each voxel based on estimated model
        gcp; H = ParforProgMon('Prewhitening: ', size(dataOut,1), 1);
        parfor pixel = 1:size(dataOut,1)    
                dataOut(pixel,:) = transpose(infer(EstMdl, double(dataOut(pixel,:))'));
                H.increment(); %#ok<PFBNS>
        end
        
        % reshape back to original data dims
        dataOut = reshape(dataOut, [n_width, n_depth, n_trials, n_windows]);        
        dataOut = permute(dataOut, [2, 1, 4, 3]); 
    end            
end