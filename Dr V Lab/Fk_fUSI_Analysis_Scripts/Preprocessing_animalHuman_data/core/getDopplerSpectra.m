% getDoppler beamforms the raw IQ data from .bin files to Doppler data (dataOut)
%
% data in:
% nChannels, int, (number of doppler channels <=> image width)
% ensLength, int, (doppler ensemble length)
% nWindows, int, (time interval to analyze [s])
% indsOfInterest, 1 x nTrials double, (fUS indices for start of each trial)
%
% data out:
% iDop{spectra} (z_pixels, x_pixels, nTimeIndices(nWindows+1), nTrials)
% spectra is a n x 2 list of filter bandwidths (in mm/s)

function [iDop, spectra] = getDopplerSpectra(nChannels, ensLength, FramerateUF, nWindows, indsOfInterest)

% error checking
if sum(isnan(indsOfInterest))~=0
    warning('NaN entries found in inds of interest!')
    indsOfInterest(isnan(indsOfInterest)) = [];
end

% velocities of interest (in mm/s)
spectra = [ 0.5 5.0; ...
    5.0 10; ...
    10 24];

% number of trials we are going to process
nTrials = length(indsOfInterest);

% set up the parallel pool
thePool = gcp; %#ok<NASGU>

% allocate memory
iDop = cell(1,size(spectra,1));

% for bandIndex = 1:size(spectra,1)
for trial = 1:nTrials
    fprintf('Beamforming: trial %i out of %i total...', trial, nTrials)
    startBin = indsOfInterest(trial); % first bin number for this trial
    for window = 1:nWindows
        % open the .bin, read it, close it
        fid = fopen(sprintf('fUS_block_%.3d.bin', startBin+window-1), 'r');
        bin_data = fread(fid, 'double');
        fclose(fid);
        
        % reshape the IQ data to 3D (2D image + complex conj)
        IQ_temp = reshape(bin_data, [],nChannels*2, ensLength);
        IQ_complex = IQ_temp(:,1:nChannels,:)+1i*IQ_temp(:,nChannels+1:2*nChannels,:);
        
        % get size of data & reshape to complex conj
        [nz, nx, nt] = size(IQ_complex);
        IQ_complex = reshape(IQ_complex, [nz*nx, nt]);
        
        % compute covariance matrix & eigenvalues
        cov_matrix = IQ_complex'*IQ_complex;
        [Eig_vect, Eig_val] = eig(cov_matrix);
        Eig_vect = fliplr(Eig_vect);
        Eig_val = rot90(Eig_val,2);
        M_A = IQ_complex*Eig_vect;
        
        % capture 85% of the variance (subtract tissue)
        skipped_eig_val = 1 : round(0.15*ensLength);
        IQF_tissu = M_A(:,skipped_eig_val)*Eig_vect(:,skipped_eig_val)';
        IQF_tissu = reshape(IQF_tissu, [nz, nx, nt]);
        IQ_complex = reshape(IQ_complex, [nz, nx, nt]);     
        IQ_complex = IQ_complex - IQF_tissu; 
        
        parfor bandIndex = 1:size(spectra,1)        
            % butterworth filter
            filterFreqs = spectra(bandIndex,:) * 20;
            [B, A] = butter(4,  filterFreqs/FramerateUF, 'bandpass'); % 20-100 is 1-5 mm/s
            IQF_corrected = filter(B, A, IQ_complex, [], 3);
            
            % compute doppler power of corrected (tissue filtered) IQ
            Dop = mean(abs(IQF_corrected).^2,3);
            iDop_tmp(bandIndex,:,:,window,trial) = interp2(Dop,2);
        end
        % this secondary assigment loop is necessary because of the way
        % that matlab parallel computing handles sliced cell array vars
        for bandIndex = 1:size(iDop_tmp,1)
            iDop{bandIndex} = squeeze(iDop_tmp(bandIndex,:,:,:,:));
        end
    end
    fprintf('done.\n')
    waitbar(trial/nTrials)
end
fprintf('Beamforming complete.\n')
end
