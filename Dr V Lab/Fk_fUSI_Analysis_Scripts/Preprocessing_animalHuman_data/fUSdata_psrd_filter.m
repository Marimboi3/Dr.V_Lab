% ***** Function corrects temporal periodic, spikes,rigid and drifts motion artifacts
% img4D = raw 4D fUS acquisition data
% params = [rgd p_rgd dtrn p_dtrn brd p_brd medn p_medn]
% p_brd = threshold number of highest frequecy peaks to remove in the
% frequency domain [ typically between 5 and 9]
% p_medn = Integer power/order to use for median smoothing filtering
%p_dtrn = Size of sliding window to use for detrending 
% isD = performs detrending first if isD = 1 
 
function imgCorrected = fUSdata_psrd_filter(img4D, params, isD)

is_rgd = params(1); par_rgd = params(2);
is_dtrn = params(3); par_dtrn = params(4);
% is_brd = params(5); par_brd = params(6);
is_lwps = params(5);  par_lwpT1 = params(6); par_lwpT2 = params(7); par_lwpFq = params(8); 
lpsParam = [par_lwpT1, par_lwpT2, par_lwpFq];
is_medn = params(9); par_medn = params(10);

disp(size(img4D));
preMC = permute(img4D, [1 3 2 4]); disp(size(preMC));


if(isD == 1)
    if(is_dtrn == 1)
        % *** Detrend correction
        disp('Start detrending drift motion....');
        fUS_dtrdC = detrend_sliding_window(preMC, par_dtrn);
        disp('Done Drift Correction !!'); disp(size(fUS_dtrdC));
        preMC = fUS_dtrdC;
    end
    if(is_rgd  == 1)
        % *** Rigid Motion Correct
        disp('Start Rigid motion correction ....');
        [rigid_fUSdataCorrected, template] = normcorre_doppler(preMC);
        disp('Rigid motion correction Done!'); disp(size(rigid_fUSdataCorrected));
        preMC = rigid_fUSdataCorrected;
    end
end

if(isD == 0)
    if(is_rgd  == 1)
        % *** Rigid Motion Correct
        disp('Start Rigid motion correction ....');
        [rigid_fUSdataCorrected, template] = normcorre_doppler(preMC);
        disp('Rigid motion correction Done!'); disp(size(rigid_fUSdataCorrected));
        preMC = rigid_fUSdataCorrected;
    end
    if(is_dtrn == 1)
        % *** Detrend correction
        disp('Start detrending drift motion....');
        fUS_dtrdC = detrend_sliding_window(preMC, par_dtrn);
        disp('Done Drift Correction !!'); disp(size(fUS_dtrdC));
        preMC = fUS_dtrdC;
    end
end

% if(is_brd == 1)
%     % *** Breathing Correction 
%     disp('Start Periodic/breathing correction ....');
%     fUS_pC =  fUS_remove_periodicMotion(preMC, par_brd);
%     disp('Done Breathing correction !!'); disp(size(fUS_pC));
%     preMC = fUS_pC;
% end

if(is_lwps == 1)
    % *** Low Pass filter to remove high frequncy noise 
    disp('Start Low pass filtreing correction ....');
    fUS_lwpC =  fUS_lowpass(preMC, lpsParam);
    disp('Done Lowpass filtering !!'); disp(size(fUS_lwpC));
    preMC = fUS_lwpC;
end

if(is_medn == 1)
    % *** Smoothing median filter 
    disp('Start Median filtering....');
    fUS_mdC = fUS_median_filter(preMC, par_medn);
    disp('Done Median filtering !!!'); disp(size(fUS_mdC));
    preMC = fUS_mdC;
end

% *** Return corrected fUS image data
imgCorrected = permute(preMC, [1 3 2 4]); disp(size(imgCorrected));
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


   function img_FC = fUS_lowpass(img4D, ftParam)
   
   preMC = permute(img4D, [1 3 2 4]); 
   disp(size(preMC));
   
    tp_img = squeeze(preMC);  
    pDat = permute(tp_img, [3 1 2]);
    
    unflt_vec = reshape(pDat, size(pDat, 1)*size(pDat, 2)*size(pDat, 3), 1);
    flt_vec = lowpass(unflt_vec, ftParam(2), ftParam(3));
%     flt_vec = bandpass(unflt_vec, [ftParam(1) ftParam(2)], ftParam(3));
    
    flt3D = reshape(flt_vec, size(pDat, 1), size(pDat, 2), size(pDat, 3));
    flt3D = permute(flt3D, [3 2 1]);
    
    tpMC = permute(flt3D, [1 4 2 3]); disp(size(tpMC));
    img_FC = permute(tpMC, [3 1 2 4]); disp(size(img_FC));
           
end
