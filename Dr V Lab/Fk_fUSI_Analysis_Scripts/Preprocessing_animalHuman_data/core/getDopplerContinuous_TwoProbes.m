% this function was adapted from code by DM
% getDoppler returns doppler signal (dataOut)
%
% data in:
% nChannels, int, (number of doppler channels <=> image width)
% ensLength, int, (doppler ensemble length)
% nBins, int, (number of .bin files/images to form)
%
% data out:
% iDop, (zPixels x xPixels(nChannels) x nTimeIndices(nBins))

function iDop = getDopplerContinuous_TwoProbes(nChannels, ensLength, FramerateUF, nBins) %#ok<INUSL>

thePool = gcp; %#ok<NASGU>
H = ParforProgMon('Beamforming Data: ', nBins, 1);
parfor i = 1:nBins
    % open the .bin, read it, close it
    fid = fopen(sprintf('fUS_block_%.3d.bin', i), 'r');
    binData = fread(fid, 'double');
    fclose(fid);
    
    % reshape the IQ data to 3D (2D image + complex conj)
    IQ_temp = reshape(binData, [], nChannels*2, ensLength);
    IQ_complex = IQ_temp(:,1:nChannels,:)+1i*IQ_temp(:,nChannels+1:2*nChannels,:);
    
    probeChannels = [1 128; 129 256];
    
    for probe = 1:2
        % Select data for each probe
        IQ_probe = IQ_complex(:,probeChannels(probe,1):probeChannels(probe,2),:);
        
        % get size of data & reshape to complex conj
        [nz, nx, nt] = size(IQ_probe);
        IQ_probe = reshape(IQ_probe, [nz*nx, nt]);
        
        % compute covariance matrix & eigenvalues
        cov_matrix = IQ_probe'*IQ_probe;
        [Eig_vect, Eig_val]= eig(cov_matrix);
        Eig_vect=fliplr(Eig_vect);
        Eig_val=rot90(Eig_val,2);
        M_A = IQ_probe*Eig_vect;
        
        % capture 85% of the variance
        skipped_eig_val = 1 : round(0.15*ensLength);
        IQF_tissu = M_A(:,skipped_eig_val)*Eig_vect(:,skipped_eig_val)';
        IQF_tissu = reshape(IQF_tissu, [nz, nx, nt]);
        IQ_probe = reshape(IQ_probe, [nz, nx, nt]);
        IQF_corrected = IQ_probe-IQF_tissu;
        
        % butterworth filter
        %             [B,A] = butter(6,  [20 100]/FramerateUF, 'bandpass'); % 20-100 is 1-5 mm/s
        %             IQF_corrected = filter(B, A, IQF_corrected, [], 3);
        %
        
        % compute doppler power of corrected (tissue filtered) IQ
        Dop = mean(abs(IQF_corrected).^2,3);
        if probe == 1
            iDop_temp = interp2(Dop,2);
        elseif probe == 2
            iDop_temp = cat(2,iDop_temp,interp2(Dop,2));
        end
    end %End for over probe
    iDop(:,:,i) = iDop_temp;
    
    % update progress bar
    H.increment();
end

end
