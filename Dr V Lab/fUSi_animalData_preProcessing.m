% ***** Function corrects temporal periodic, spikes,rigid and drifts motion artifacts
% img4D = raw 4D fUS acquisition data
% params = [rgd p_rgd dtrn p_dtrn brd p_brd medn p_medn]
 
function imgCorrected = fUSi_animalData_preProcessing(img4D, params, bsTempt)

is_rgd = params{1}{1}; 
refLen = params{1}{2}; 
binWidth = params{1}{3};

is_lwps = params{2}{1};
par_lwpT1 = params{2}{2}; 
par_lwpT2 = params{2}{3}; 
par_lwpFq = params{2}{4}; 
lpsParam = [par_lwpT1, par_lwpT2, par_lwpFq];

disp(size(img4D));
preMC = permute(img4D, [1 3 2 4]); disp(size(preMC));

    if(is_rgd  == 1)
        % *** Rigid Motion Correct
        disp('Start Rigid motion correction ....');
        
        if ~exist('bsTempt','var')
            [rigid_fUSdataCorrected, template] = normcorre_doppler_human_preProcessing_ka(preMC, refLen, binWidth);
        else
            baselineTemplate = bsTempt;
             [rigid_fUSdataCorrected, template] = normcorre_doppler_human_preProcessing_ka(preMC, refLen, binWidth, baselineTemplate);
        end
               
        disp('Rigid motion correction Done!'); disp(size(rigid_fUSdataCorrected));
        preMC = rigid_fUSdataCorrected;
    end
    
    if(is_lwps == 1)
        % *** Low Pass filter to remove high frequncy noise 
        disp('Start Low pass filtreing correction ....');
        fUS_lwpC =  fUS_lowpass(preMC, lpsParam);
        disp('Done Lowpass filtering !!'); disp(size(fUS_lwpC));
        preMC = fUS_lwpC;
    end

% *** Return corrected fUS image data
imgCorrected = permute(preMC, [1 3 2 4]); disp(size(imgCorrected));
disp('Done fUS data periodic, rigid, smooting and drift motion correction !!!!!');

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

