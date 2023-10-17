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

function iDop = getDopplerContinuous(nChannels, ensLength, FramerateUF, nBins) %#ok<INUSL>

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
    
    % get size of data & reshape to complex conj
    [nz, nx, nt] = size(IQ_complex);
    IQ_complex = reshape(IQ_complex, [nz*nx, nt]);
    
    % compute covariance matrix & eigenvalues
    cov_matrix = IQ_complex'*IQ_complex;
    [Eig_vect, Eig_val] = eig(cov_matrix);
    Eig_vect = fliplr(Eig_vect);
    Eig_val = rot90(Eig_val,2); %#ok<NASGU>
    M_A = IQ_complex*Eig_vect;
    
    % capture 85% of the variance (?)
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
    iDop(:,:,i) = mean(abs(IQF_corrected).^2,3);
    
    % interpolate image for smoothness
    %             iDop(:,:,i) = interp2(Dop,2);    
    % note: I am removing this Oct 29, 2020 to reduce file size.
    
    % update progress bar
    H.increment();
end

end
