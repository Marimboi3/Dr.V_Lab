% test script

%% hard coded params
nChannels = 128;
ensLength = 500;
FramerateUF = 500;
nWindows = 10;
indsOfInterest = [518 528 550 600];

cd 'Z:\fUS\Raw Data\20170921\Loop_fUS_2017-09-21@17-13-57'

%% testing serial

% number of trials we are going to process
nTrials = length(indsOfInterest);

for trial = 1:nTrials
    startBin = indsOfInterest(trial); % first bin number for this trial
    fprintf('Beamforming: %i trial out of %i total...', trial, nTrials)
    
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
        [B,A] = butter(4,  [20 100]/FramerateUF, 'bandpass'); % 20-100 is 1-5 mm/s
        IQF_corrected = filter(B, A, IQF_corrected, [], 3);
        
        % compute doppler power of corrected (tissue filtered) IQ
        Dop = mean(abs(IQF_corrected).^2,3);
        iDop(:,:,window,trial) = interp2(Dop,2);
    end
    fprintf('done.\n')
end
fprintf('Beamforming complete.\n')

%% testing GPU
clear iDop

% number of trials we are going to process
nTrials = length(indsOfInterest);

for trial = 1:nTrials
    startBin = indsOfInterest(trial); % first bin number for this trial
    fprintf('Beamforming: %i trial out of %i total...', trial, nTrials)
    
    for window = 1:nWindows
        % open the .bin, read it, close it
        fid = fopen(sprintf('fUS_block_%.3d.bin', startBin+window-1), 'r');
        bin_data = fread(fid, 'double');
        fclose(fid);
        
        % convert to GPU array
        bin_data = gpuArray(bin_data);
        
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
        [B, A] = butter(4,  [20 100]/FramerateUF, 'bandpass'); % 20-100 is 1-5 mm/s
        IQF_corrected = filter(B, A, IQF_corrected, [], 3);
                
        % compute doppler power of corrected (tissue filtered) IQ
        Dop = mean(abs(IQF_corrected).^2,3);
        iDop(:,:,window,trial) = interp2(Dop,2);
    end
    fprintf('done.\n')
end
fprintf('Beamforming complete.\n')

%% testing parallel

% number of trials we are going to process
nTrials = length(indsOfInterest);

% set up the parallel pool
thePool = gcp;

for trial = 1:nTrials
    startBin = indsOfInterest(trial); % first bin number for this trial
    fprintf('Beamforming: %i trial out of %i total...', trial, nTrials)
    
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
        [B,A] = butter(4,  [20 100]/FramerateUF, 'bandpass'); % 20-100 is 1-5 mm/s
        IQF_corrected = filter(B, A, IQF_corrected, [], 3);
        
        
        % compute doppler power of corrected (tissue filtered) IQ
        Dop = mean(abs(IQF_corrected).^2,3);
        iDop(:,:,window,trial) = interp2(Dop,2);
    end
    fprintf('done.\n')
end
fprintf('Beamforming complete.\n')


