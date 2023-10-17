%% Data Analysis: 20201210 ICA Oclusion ByPass PN stimulation

% *** Init 
clear all; close all; 
clearvars; clc;

bp_acqPeriod = 0.4;
bp_acqFreq = 1 / bp_acqPeriod;

% *** Load Motion Corrected Data ******
Fn_byps_mc = 'bypass12_10_mc.acq';
load_byps = load(Fn_byps_mc, '-mat');
bypsFull_Data = squeeze(load_byps.Acquisition.Data);   
bypsFull_times = load_byps.Acquisition.T;

byps_FullData_avg = squeeze(mean(permute(bypsFull_Data, [2 1 3]), 3));

%  *** Display Average Image ***
figure; 
colormap(turbo); imagesc(byps_FullData_avg); title('Right MVD PNS fUS Sequence Average');

% *** ROI averaged data *****
rois_byps = '20201210_Bypass12_10_mc_1_11.txt';
load_byps = importdata(rois_byps);

byps_roi_times = load_byps.data(:, 1);
byps_roi_1 = load_byps.data(:, 2);
byps_roi_2 = load_byps.data(:, 3);
byps_roi_3 = load_byps.data(:, 4);
byps_roi_4 = load_byps.data(:, 5);
byps_roi_5 = load_byps.data(:, 6);
byps_roi_6 = load_byps.data(:, 7);
byps_roi_7 = load_byps.data(:, 8);
byps_roi_8 = load_byps.data(:, 9);
byps_roi_9 = load_byps.data(:, 10);
byps_roi_10 = load_byps.data(:, 11);
byps_roi_11 = load_byps.data(:, 12);
% byps_roi_12 = load_byps.data(:, 13);
% byps_roi_13 = load_byps.data(:, 14);
% byps_roi_14 = load_byps.data(:, 15);
% byps_roi_15 = load_byps.data(:, 16);

nRois = 11;

  % *** Display FUS ROIs time Data and Stim ON/OFF Results *****
for nr = 1 : nRois
    
    % ***** Trial Stimulation ON and OFF Epoch Matices *****
    byps_roi_stim_data = eval(strcat('byps_roi_', num2str(nr)));
    byps_roi_title =  strcat('ByPass PN Stimulation on/off Activation (Exp 6):  ROI- ', num2str(nr));
    byps_roi_hdr = strcat('Temporal Plot:  ROI- ', num2str(nr));

    bp_roi_bsline = mean(byps_roi_stim_data(100 : 160));
%     byps_roi_stim_data_flt = byps_roi_stim_data;
%     byps_roi_stim_data_flt(1) = mean(byps_roi_stim_data(3:8));
%     byps_roi_stim_data_flt(length(byps_roi_stim_data)) = mean(byps_roi_stim_data(length(byps_roi_stim_data)-8 : length(byps_roi_stim_data)) - 3);
    
    byps_roi_stim_data_norm = 100*((byps_roi_stim_data - bp_roi_bsline) / bp_roi_bsline);
    
    byps_roi_bsLine = (byps_roi_stim_data(131 : 230) - bp_roi_bsline) / bp_roi_bsline;
    byps_roi_stimON = 100*((byps_roi_stim_data(231 : 503) - bp_roi_bsline) / bp_roi_bsline);
    byps_roi_stimOFF = 100*((byps_roi_stim_data(504 : 595) - bp_roi_bsline) / bp_roi_bsline);
    
        bp_off = 1 : 92; bp_on = 93:365;
        bp_roi_tm = (1 : size(byps_roi_stim_data, 1)) * bp_acqPeriod;
        
        % *** Display Results
        figure; hold on;
        plot(bp_roi_tm, byps_roi_stim_data_norm); title(byps_roi_hdr); ylabel('PD (a.u.)'); xlabel('time (s)');
        xline(1,'--b',{'PreStim'}); xline(160,'--g',{'StimON'}); xline(230,'--k',{'StimEff'}); xline(503,'--r',{'StimOFF'}); ylabel('% Change');
        hold off;
        
%         figure; hold on;
%         subplot(212); plot(byps_roi_stim_data); title(byps_roi_hdr); ylabel('PD (a.u.)'); xlabel('time (s)');
%         xline(1,'--b',{'PreStim'}); xline(160,'--g',{'StimON'}); xline(230,'--k',{'StimEff'}); xline(503,'--r',{'StimOFF'}); ylabel('PD (a.u.) ');
%         subplot(211);
%         shadedErrorBar(bp_off, byps_roi_stimOFF, {@mean,@std},'lineprops','-b','transparent',1, 'patchSaturation',0.1);
%         shadedErrorBar(bp_on, byps_roi_stimON, {@mean,@std},'lineprops','-r','transparent',1, 'patchSaturation',0.1);
%         xline(1,'--r',{'StimOFF'}); 
%         xline(92,'--g',{'StimON'});
%         title(strcat('  ', byps_roi_title));
%         xlabel('time (s)'); ylabel('% Chanage');
%         hold off;
end


% ****  T-test to determine pixels of significant differences *****
byps_bsline = squeeze(mean(mean(mean(bypsFull_Data(:, :, 190 : 230)))));

byps_preStim = squeeze(permute(bypsFull_Data(:, :, (131:230)), [2 1 3]));
byps_stimON = squeeze(permute(bypsFull_Data(:, :, (231 : 503)), [2 1 3]));
byps_stimON_early = squeeze(permute(bypsFull_Data(:, :, (231:330)), [2 1 3]));
byps_stimON_late = squeeze(permute(bypsFull_Data(:, :, (404:503)), [2 1 3]));
byps_stimOFF = squeeze(permute(bypsFull_Data(:, :, (509 : 600)), [2 1 3]));

byps_data_avg = squeeze(mean(permute(bypsFull_Data, [2 1 3]), 3));

% **** 
[eStOn_rject, eStOn_p_vals] = ttest2(permute(byps_preStim, [3 1 2]), permute(byps_stimON_early, [3 1 2]));
[lStOn_rject, lStOn_p_vals] = ttest2(permute(byps_preStim, [3 1 2]), permute(byps_stimON_late, [3 1 2]));
[StOn_rject, StOn_p_vals] = ttest2(permute(byps_preStim, [3 1 2]), permute(byps_stimON, [3 1 2]));
[StOff_rject, StOff_p_vals] = ttest2(permute(byps_preStim, [3 1 2]), permute(byps_stimOFF, [3 1 2]));

eStOn_p_vals  = squeeze(eStOn_p_vals);
lStOn_p_vals  = squeeze(lStOn_p_vals);
StOn_p_vals  = squeeze(StOn_p_vals);
StOff_p_vals  = squeeze(StOff_p_vals);

bp_sigAlpha = 1e-10;
eStOn_sig_indx = -log(eStOn_p_vals .* (eStOn_p_vals < bp_sigAlpha));
LStOn_sig_indx = -log(lStOn_p_vals .* (lStOn_p_vals < bp_sigAlpha));
StOn_sig_indx = -log(StOn_p_vals .* (StOn_p_vals < bp_sigAlpha));

% eStOn_sig_indx(find(eStOn_sig_indx == inf)) = 0.0;
LStOn_sig_indx(find(LStOn_sig_indx== inf)) = 0.0;
% StOn_sig_indx(find(StOn_sig_indx == inf)) = 0.0;

% **** Significant ROIs
% bp_sig_roi_9 = roipoly(LStOn_sig_indx);
% bpRoi_all = 1.0*bp_sig_roi_1 + 1.0*bp_sig_roi_2 + 1.0*bp_sig_roi_3 + 1.0*bp_sig_roi_4 + 1.0*bp_sig_roi_5 + 1.0*bp_sig_roi_6 + 1.0*bp_sig_roi_7 + 1.0*bp_sig_roi_8 + 1.0*bp_sig_roi_9;
% bpSig_roi_1 = squeeze(mean(mean(1.0*bp_sig_roi_1 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_2 = squeeze(mean(mean(1.0*bp_sig_roi_2 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_3 = squeeze(mean(mean(1.0*bp_sig_roi_3 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_4 = squeeze(mean(mean(1.0*bp_sig_roi_4 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_5 = squeeze(mean(mean(1.0*bp_sig_roi_5 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_6 = squeeze(mean(mean(1.0*bp_sig_roi_6 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_7 = squeeze(mean(mean(1.0*bp_sig_roi_7 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_8 = squeeze(mean(mean(1.0*bp_sig_roi_8 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_9 = squeeze(mean(mean(1.0*bp_sig_roi_9 .* permute(bypsFull_Data, [2 1 3]))));
% bpSig_roi_10 = squeeze(mean(mean(1.0*bp_sig_roi_10 .* permute(bypsFull_Data, [2 1 3]))));

