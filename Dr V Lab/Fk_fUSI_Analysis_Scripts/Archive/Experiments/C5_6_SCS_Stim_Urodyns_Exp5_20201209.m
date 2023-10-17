%% Data Analysis: 20201209 C5-6 Spinal Cord Imaging and Urodynamics
% No stimulation

% *** Init 
% clear all; close all; 
% clearvars; clc;

% *** Load Motion Corrected Data ******
Fn_uds_raw = 'UDS filling.acq';
Fn_uds_mc = 'UDS_filling_mc_1_21_0.acq';

% *** ROIS ****
load('20201209_UDC_rois_sigPvals_selected_ROIs_logic.mat');

load_uds = load(Fn_uds_mc, '-mat');
uds_Data = squeeze(load_uds.Acquisition.Data);   
uds_times = load_uds.Acquisition.T;

uds_bsline = squeeze(mean(mean(mean(uds_Data(:, :, 59 : 119)))));

uds_empty = squeeze(permute(uds_Data(:, :, (1:120)), [2 1 3]));
uds_flng_early = squeeze(permute(uds_Data(:, :, (121:240)), [2 1 3]));
uds_flng_late = squeeze(permute(uds_Data(:, :, (361:480)), [2 1 3]));
uds_full = squeeze(permute(uds_Data(:, :, (481 : 600)), [2 1 3]));

uds_data_avg = squeeze(mean(permute(uds_Data, [2 1 3]), 3));

% **** T-test to determine pixels of significant differences *****
[eF_rject, emptyFull_p_vals] = ttest2(permute(uds_empty, [3 1 2]), permute(uds_full, [3 1 2]));
[eEF_rject, emptyEly_p_vals] = ttest2(permute(uds_empty, [3 1 2]), permute(uds_flng_early, [3 1 2]));
[eLF_rject, emptyLate_p_vals] = ttest2(permute(uds_empty, [3 1 2]), permute(uds_flng_late, [3 1 2]));

emptyFull_p_vals = squeeze(emptyFull_p_vals);
emptyEly_p_vals = squeeze(emptyEly_p_vals);
emptyLate_p_vals  = squeeze(emptyLate_p_vals);

sigAlpha = 1e-25;
emptyFull_sig_indx = -log(emptyFull_p_vals .* (emptyFull_p_vals < sigAlpha));
emptyEly_sig_indx = -log(emptyEly_p_vals .* (emptyEly_p_vals < sigAlpha));
emptyLate_sig_indx = -log(emptyLate_p_vals .* (emptyLate_p_vals < sigAlpha));

% emptyFull_sig_indx(find(emptyFull_sig_indx == inf)) = 0.0;
% emptyEly_sig_indx(find(emptyEly_sig_indx == inf)) = 0.0;
%  emptyLate_sig_indx(find(emptyLate_sig_indx == inf)) = 0.0;

% % **** Extract Significant ROIs *****
% uds_sig_roi_20 = roipoly(emptyFull_sig_indx);
% udsRoi_all = 1.0*uds_sig_roi_1 + 1.0*uds_sig_roi_2 + 1.0*uds_sig_roi_3 + 1.0*uds_sig_roi_4 + 1.0*uds_sig_roi_5 + 1.0*uds_sig_roi_6 
% + 1.0*uds_sig_roi_7 + 1.0*uds_sig_roi_8 + 1.0*uds_sig_roi_9 + 1.0*uds_sig_roi_10 + 1.0*uds_sig_roi_11 + 1.0*uds_sig_roi_12 + 1.0*uds_sig_roi_13; 
udsSig_roi_1 = squeeze(mean(mean(uds_sig_roi_1 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_2 = squeeze(mean(mean(uds_sig_roi_2 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_3 = squeeze(mean(mean(uds_sig_roi_3 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_4 = squeeze(mean(mean(uds_sig_roi_4 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_5 = squeeze(mean(mean(uds_sig_roi_5 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_6 = squeeze(mean(mean(uds_sig_roi_6 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_7 = squeeze(mean(mean(uds_sig_roi_7 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_8 = squeeze(mean(mean(uds_sig_roi_8 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_9 = squeeze(mean(mean(uds_sig_roi_9 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_10 = squeeze(mean(mean(uds_sig_roi_10 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_11 = squeeze(mean(mean(uds_sig_roi_11 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_12 = squeeze(mean(mean(uds_sig_roi_12 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_13 = squeeze(mean(mean(uds_sig_roi_13 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_14 = squeeze(mean(mean(uds_sig_roi_14 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_15 = squeeze(mean(mean(uds_sig_roi_15 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_16 = squeeze(mean(mean(uds_sig_roi_16 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_17 = squeeze(mean(mean(uds_sig_roi_17 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_18 = squeeze(mean(mean(uds_sig_roi_18 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_19 = squeeze(mean(mean(uds_sig_roi_19 .* permute(uds_Data, [2 1 3]))));
udsSig_roi_20 = squeeze(mean(mean(uds_sig_roi_20 .* permute(uds_Data, [2 1 3]))));

nUds_rois = 20;
udsRoi_all = 1.0*uds_sig_roi_1;

for nrs = 2 : nUds_rois
    tp_roi = eval(strcat('uds_sig_roi_', num2str(nrs)));
    udsRoi_all = udsRoi_all + 1.0*tp_roi;
end

% *** Display  Pre Resection Results ****
   for nuds = 1 : nUds_rois
       
        Uds_title =  strcat(' Urodynamics fUS Regions with p-vals ROI- ', num2str(nuds));
        Uds_stim_data = eval(strcat('udsSig_roi_', num2str(nuds)));
        Uds_bsline = mean(Uds_stim_data(160 : 180));
        
        Uds_stim_data_norm = 100*((Uds_stim_data - Uds_bsline) / Uds_bsline);
        
%         % *** Display Results
%         figure; 
%         hold on; 
%         plot(uds_times, Uds_stim_data);
%         xlabel('time (s)'); ylabel('PD a.u.'); title(Uds_title); xline(1,'--k',{'empty'}); xline(180,'--g',{'filling'}); xline(540,'--r',{'Full'});
%         hold off;
        
        % *** Display as percentage change
        figure; 
        hold on; 
        plot(uds_times, Uds_stim_data_norm);
        xlabel('time (s)'); ylabel('% change'); title(Uds_title); xline(1,'--k',{'empty'}); xline(180,'--g',{'filling'}); xline(540,'--r',{'Full'});
        hold off;
        
   end
   
   
   %  *** Display Average Image ***
figure; 
subplot(121); colormap(turbo); imagesc(uds_data_avg); title('C5-6 Spinal Cord Imaging and Urodynamics fUS Sequence Average');
% subplot(122); colormap(hot); imagesc(udsRoi_all); title('Selected ROIs with Significant p-vals');
   
figure; % Average fUS Image and Significant p-vals
subplot(131);  colormap(hot); imagesc(emptyEly_sig_indx); title(strcat('Empty and Early Filling - (-log(p - vals < ', string(sigAlpha), '))')); colorbar;
subplot(132); imagesc(emptyLate_sig_indx); title(strcat('Empty and Late Filling - (-log(p - vals < ', string(sigAlpha), '))')); colorbar;
subplot(133); imagesc(emptyFull_sig_indx); title(strcat('Empty and Full - (-log(p - vals < ', string(sigAlpha), '))')); colorbar;

% figure;
% imagesc(emptyFull_sig_indx); title(strcat('Empty and Full - (-log(p - vals < ', string(sigAlpha), '))')); colorbar;
% 
figure; % Selected ROIs
colormap(gray); imagesc(udsRoi_all); title('Selected ROIs with Significant p-vals');
