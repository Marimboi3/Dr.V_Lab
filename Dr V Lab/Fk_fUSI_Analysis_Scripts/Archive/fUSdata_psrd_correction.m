% ***** Function corrects temporal periodic, spikes,rigid and drifts motion artifacts
% img4D = raw 4D fUS acquisition data
% nPks = threshold number of highest frequecy peaks to remove in the
% frequency domain [ typically between 5 and 9]
% nOrder = Integer power/order to use for median smoothing filtering
% sWdow = Size of sliding window to use for detrending 
% idD = performs detrending first if isD = 1 
function imgCorrected = fUSdata_psrd_correction(img4D, nPks, nOrder, sWdow, isD)
disp(size(img4D));
preMC = permute(img4D, [1 3 2 4]); disp(size(preMC));


if(isD == 1)
    % *** Detrend correction
    disp('Start detrending drift motion....');
    fUS_dtrdC = detrend_sliding_window(preMC, sWdow);
    disp('Done Drift Correction !!'); disp(size(fUS_dtrdC));
    preMC = fUS_dtrdC;
end

% *** Rigid Motion Correct
disp('Start Rigid motion correction ....');
[rigid_fUSdataCorrected, template] = normcorre_doppler(preMC);
disp('Rigid motion correction Done!'); disp(size(rigid_fUSdataCorrected));

if(isD == 0)
    % *** Detrend correction
    disp('Start detrending drift motion....');
    fUS_dtrdC = detrend_sliding_window(rigid_fUSdataCorrected, sWdow);
    disp('Done Drift Correction !!'); disp(size(fUS_dtrdC));
    rigid_fUSdataCorrected = fUS_dtrdC;
end

% *** Breathing Corrention 
disp('Start Periodic/breathing correction ....');
fUS_pC =  fUS_remove_periodicMotion(rigid_fUSdataCorrected, nPks);
disp('Done Breathing correction !!'); disp(size(fUS_pC));

% *** Smoothing median filter 
disp('Start Median filtering....');
fUS_mdC = fUS_median_filter(fUS_pC, nOrder);
disp('Done Median filtering !!!'); disp(size(fUS_mdC));

% *** Return corrected fUS image data
imgCorrected = permute(fUS_mdC, [1 3 2 4]); disp(size(imgCorrected));
disp('Done fUS data periodic, rigid, smooting and drift motion correction !!!!!');

end


% **** Given temporal 4D fUS data and number of frequency peaks, removes top k
% **** frequencies in the freq domain
% img4D = 4D fUS acquisition data
% nPks = threshold number of highest frequecy peaks to remove in the
function img_FC = fUS_remove_periodicMotion(img4D, fct)

    tp_img = squeeze(img4D);
    fUS_FC = [];
    for yj = 1 : size(tp_img, 1)
        for xi = 1 : size(tp_img, 2)
            nxt_px = squeeze(tp_img(yj, xi, :));
            nxt_px_fft = fftshift(fft(nxt_px));

            [top_mx, top_indx] = maxk(nxt_px_fft, fct);
            disp(top_indx);
            if(top_indx > 1 & top_indx < size(nxt_px_fft,1))
                    nxt_px_fft(top_indx) = (nxt_px_fft(top_indx - 1) + nxt_px_fft(top_indx + 1))/2;
            end

%             top_indx_wo_mx = top_indx(2 : length(top_indx));
%             nxt_px_fft(top_indx) = 0;

            nxt_inpx = abs(ifft(nxt_px_fft));
            fUS_FC(yj, xi, :) = nxt_inpx;
        end
    end
    tp_fn = permute(fUS_FC, [1 4 2 3]);
    img_FC = permute(tp_fn, [1 3 2 4]);
%     img_FC = permute(fUS_FC, [1 4 2 3]);
    
end

% **** Given temporal 3D fU data and number/order, applies a median filter
% **** with order given in time domain
% img4D = 4D fUS acquisition data
% nPks = Integer power/order to use for median smoothing filtering
function img_FC = fUS_median_filter(img4D, odr_k)

    tp_img = squeeze(img4D);
    fUS_FC = [];
    
    for yj = 1 : size(tp_img, 1)
        for xi = 1 : size(tp_img, 2)
            nxt_px = squeeze(tp_img(yj, xi, :));
            nxt_px_md = medfilt1(nxt_px, odr_k);
            fUS_FC(yj, xi, :) = nxt_px_md;
        end
    end
    tp_fn = permute(fUS_FC, [1 4 2 3]);
%     img_FC = permute(fUS_FC, [1 4 2 3]);
     img_FC = permute(tp_fn, [1 3 2 4]);
    
end
