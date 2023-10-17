% getDoppler beamforms the raw IQ data from .bin files to Doppler data (dataOut)
%
% data in:
% nChannels, int, (number of doppler channels <=> image width)
% ensLength, int, (doppler ensemble length)
% nWindows, int, (time interval to analyze [s])
% indsOfInterest, 1 x nTrials double, (fUS indices for start of each trial)
%
% data out:
% iDop, (z_pixels, x_pixels, nTimeIndices(nWindows+1), nTrials)

function iDop = getDoppler(nChannels, ensLength, FramerateUF, nWindows, indsOfInterest)

% In case nWindows is not an int, make it an int.
if ~isa(nWindows,'integer')
   warning('nWindows is not an int. Rounding number down to closet int')
   nWindows = floor(nWindows);
end

if sum(isnan(indsOfInterest))~=0
    warning('NaN entries found in inds of interest!')
    indsOfInterest(isnan(indsOfInterest)) = [];
end

% number of trials we are going to process
nTrials = length(indsOfInterest);

% set up the parallel pool
thePool = gcp; %#ok<NASGU>

% rudimentary (but robust) progress bar
fprintf(['Beamforming Progress: ' repmat('|',[1 100]) '\nCurrent Progress:     '])
progUpdate = ceil(linspace(1,nTrials));

for trial = 1:nTrials    
    startBin = indsOfInterest(trial); % first bin number for this trial        
    if sum(progUpdate==trial); fprintf('|'); end % update progress bar
    
    parfor window = 1:nWindows
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
        [Eig_vect, Eig_val]= eig(cov_matrix);
        Eig_vect=fliplr(Eig_vect);
        Eig_val=rot90(Eig_val,2);
        M_A = IQ_complex*Eig_vect;
        
        % capture 85% of the variance
        skipped_eig_val = 1 : round(0.15*ensLength);
        IQF_tissu = M_A(:,skipped_eig_val)*Eig_vect(:,skipped_eig_val)';
        IQF_tissu = reshape(IQF_tissu, [nz, nx, nt]);
        IQ_complex = reshape(IQ_complex, [nz, nx, nt]);
        IQF_corrected = IQ_complex-IQF_tissu;
        
        % butterworth filter
        %             [B,A] = butter(6,  [20 100]/FramerateUF, 'bandpass'); % 20-100 is 1-5 mm/s
        %             IQF_corrected = filter(B, A, IQF_corrected, [], 3);
        %
        
        % compute doppler power of corrected (tissue filtered) IQ
        iDop(:,:,window,trial) = mean(abs(IQF_corrected).^2,3);   
        
        % interpolate image for smoothness
        %         iDop(:,:,window,trial) = interp2(Dop,2);                
        % note: I am removing this Oct 29, 2020 to reduce file size.
    end        
end

fprintf('\nBeamforming complete.\n')
end
