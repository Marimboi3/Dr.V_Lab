%% Data Analysis - 20201118 PNS and SCS Spinal cord with equal 30s and 60s trial epoch (stim on + off) length

% *** Init 
clear all; close all; 
clearvars; clc;

% ********* Acquisition Parameters **********
acqPeriod = 1;                                                       % fUS acquisition period (s)
acqFreq = 1 / acqPeriod;                                         % Hz 

% *** Load Motion Corrected Data ******
Fn_pns_mc = '20201118_PNS_Stimulation_10_minutes_mc_1_11_0.acq';
load_pns = load(Fn_pns_mc, '-mat');
pnsFull_Data = squeeze(load_pns.Acquisition.Data);   
pnsFull_times = load_pns.Acquisition.T;

Fn_scs_mc = '20201118_SCS_Stimulation_10_minutes_mc_1_11_0.acq';
load_scs = load(Fn_scs_mc, '-mat');
scsFull_Data = squeeze(load_scs.Acquisition.Data);   
scsFull_times = load_scs.Acquisition.T;

pns_FullData_avg = squeeze(mean(permute(pnsFull_Data, [2 1 3]), 3));
scs_FullData_avg = squeeze(mean(permute(scsFull_Data, [2 1 3]), 3));

%  *** Display Average Image ***
figure; 
subplot(121);
colormap(turbo); imagesc(pns_FullData_avg); title('GekoStim PNS fUS Average');
subplot(122); imagesc(scs_FullData_avg); title('GekoStim SCS fUS Average');

% % **** Load ROI of Motion Corrected fUS Data
% rois_pns = '20201118_PNS_Stim_ROI_Data_1_11_0_v3.txt';
% rois_scs = '20201118_SCS_Stim_ROI_Data_1_11_0_v3.txt';
% load_pns = importdata(rois_pns);
% load_scs = importdata(rois_scs);
% 
% % *** PNS rois temporal data ****
% pns_roi_times = load_pns.data(:, 1);
% pns_roi_1 = load_pns.data(:, 2);
% pns_roi_2 = load_pns.data(:, 3);
% pns_roi_3 = load_pns.data(:, 4);
% pns_roi_4 = load_pns.data(:, 5);
% pns_roi_5 = load_pns.data(:, 6);
% pns_roi_6 = load_pns.data(:, 7);
% pns_roi_7 = load_pns.data(:, 8);
% pns_roi_8 = load_pns.data(:, 9);
% pns_roi_9 = load_pns.data(:, 10);
% pns_roi_10 = load_pns.data(:, 11);
% pns_roi_11 = load_pns.data(:, 12);
% pns_roi_12 = load_pns.data(:, 13);
% pns_roi_13 = load_pns.data(:, 14);
% pns_roi_14 = load_pns.data(:, 15);
% pns_roi_15 = load_pns.data(:, 16);
% % pns_roi_16 = load_pns.data(:, 17);
% % pns_roi_17 = load_pns.data(:, 18);
% % pns_roi_18 = load_pns.data(:, 19);
% % pns_roi_19 = load_pns.data(:, 20);
% % pns_roi_20 = load_pns.data(:, 21);
% % pns_roi_21 = load_pns.data(:, 22);
% % pns_roi_22 = load_pns.data(:, 23);
% % pns_roi_23 = load_pns.data(:, 24);
% % pns_roi_24 = load_pns.data(:, 25);
% nPns = 15;
% 
% % *** SCS rois temporal data****
% scs_roi_times = load_scs.data(:, 1);
% scs_roi_1 = load_scs.data(:, 2);
% scs_roi_2 = load_scs.data(:, 3);
% scs_roi_3 = load_scs.data(:, 4);
% scs_roi_4 = load_scs.data(:, 5);
% scs_roi_5 = load_scs.data(:, 6);
% scs_roi_6 = load_scs.data(:, 7);
% scs_roi_7 = load_scs.data(:, 8);
% scs_roi_8 = load_scs.data(:, 9);
% scs_roi_9 = load_scs.data(:, 10);
% scs_roi_10 = load_scs.data(:, 11);
% scs_roi_11 = load_scs.data(:, 12);
% scs_roi_12 = load_scs.data(:, 13);
% scs_roi_13 = load_scs.data(:, 14);
% scs_roi_14 = load_scs.data(:, 15);
% scs_roi_15 = load_scs.data(:, 16);
% scs_roi_16 = load_scs.data(:, 17);
% scs_roi_17 = load_scs.data(:, 18);
% scs_roi_18 = load_scs.data(:, 19);
% scs_roi_19 = load_scs.data(:, 20);
% scs_roi_20 = load_scs.data(:, 21);
% scs_roi_21 = load_scs.data(:, 22);
% scs_roi_22 = load_scs.data(:, 23);
% scs_roi_23 = load_scs.data(:, 24);
% scs_roi_24 = load_scs.data(:, 25);
% scs_roi_25 = load_scs.data(:, 26);
% scs_roi_26 = load_scs.data(:, 27);
% scs_roi_27= load_scs.data(:, 28);
% scs_roi_28 = load_scs.data(:, 29);
% scs_roi_29 = load_scs.data(:, 30);
% scs_roi_30 = load_scs.data(:, 31);
% scs_roi_31 = load_scs.data(:, 32);
% scs_roi_32 = load_scs.data(:, 33);
% 
% nScs = 32;
% 
% % ***** Generate SCS Stimulation ROI Result Plots ******8
% % **** Stimulation Parameters ******
% len_StimON_pns = 60;
% len_StimOFF_pns = 60;
% pns_section_len = len_StimON_pns + len_StimOFF_pns;                           % length of Stimulation ON and OFF
% 
% % ***** Trial Stimulation ON and OFF Epoch Matices *****
% roi_time = pns_roi_times;
% stim_section_times = roi_time; 
% stim_sections = 1 + floor((stim_section_times)/pns_section_len);
% nTrials = stim_sections(end);
% 
% ONnTpts = len_StimON_pns * acqFreq;
% OFFnTpts = len_StimOFF_pns * acqFreq;
% 
%  pns_stim_offset = 0;
% 
% for np = 1 : nPns
%         PNS_title =  strcat('20201118 PNS Stim on/off Activation - ROI-', num2str(np));
%         PNS_roi_hdr  = strcat('Temporal Plot:  ROI- ', num2str(np));
%         PNS_stim_data  = eval(strcat('pns_roi_', num2str(np)));
%         
%         pns_shft_stimData = fUS_data_shifted(PNS_stim_data, pns_stim_offset)';
%         [pns_stimON_trials_avg, pns_stimOFF_trials_avg] = get_stimONOff_trials_vec(pns_shft_stimData, stim_section_times, ONnTpts, OFFnTpts, pns_section_len);
% 
%         pns_bsline_ref = mean(mean(pns_stimOFF_trials_avg(50:59, :)));
%         pns_stimON_trials_norm = 100*(pns_stimON_trials_avg - pns_bsline_ref)/pns_bsline_ref;
%         pns_stimOFF_trials_norm = 100*(pns_stimOFF_trials_avg - pns_bsline_ref)/pns_bsline_ref;
%         
%         % *** Display Results
%         t_off_pns = 1 : size(pns_stimON_trials_avg, 1);
%         t_on_pns = size(pns_stimON_trials_avg, 1) : size(pns_stimON_trials_avg, 1) + size(pns_stimOFF_trials_avg, 1) - 1;
% 
%         figure; 
%         subplot(212); plot(PNS_stim_data); title(PNS_roi_hdr); ylabel('PD (a.u.)'); xlabel('time (s)');
%         for sm = 1 : nTrials * 2
%              if(mod(sm, 2) == 1)
%                 xoff = (sm - 1)*len_StimOFF_pns + (-1)*pns_stim_offset;
%                 if(xoff >= 0)
%                     xline(xoff,'--r',{'Off'}); 
%                 end        
%              else
%                  xon = (sm - 1)*len_StimON_pns + (-1)*pns_stim_offset;
%             if(xon >= 0)
%                 xline(xon,'--g',{'ON'});
%             end
%              end
%         end
%         
%         subplot(211); hold on;
%         shadedErrorBar(t_off_pns, pns_stimOFF_trials_norm', {@mean,@std},'lineprops','-b','transparent',1, 'patchSaturation',0.1);
%         shadedErrorBar(t_on_pns, pns_stimON_trials_norm', {@mean,@std},'lineprops','-r','transparent',1, 'patchSaturation',0.1);
%         xline(1,'--k',{'StimOFF'}); 
%         xline(len_StimON_pns,'--k',{'StimON'});
%         title(strcat('  ', PNS_title));
%         xlabel('time (s)');
%         ylabel('% Change');
%         hold off;
%         
% end
% 
% 
% % **** Generate SCS Stimulation ROI Result Plots ******
% len_StimON_scs = 30;
% len_StimOFF_scs = 30;
% scs_section_len = len_StimON_scs + len_StimOFF_scs;                           % length of Stimulation ON and OFF
% 
% % ***** Trial Stimulation ON and OFF Epoch Matices *****
% roi_time = scs_roi_times;
% stim_section_times = roi_time; 
% stim_sections = 1 + floor((stim_section_times)/scs_section_len);
% nTrials = stim_sections(end);
% 
% ONnTpts = len_StimON_scs * acqFreq;
% OFFnTpts = len_StimOFF_scs * acqFreq;
% 
%  scs_stim_offset = -20;
% 
% for ns = 1 : nScs
%         SCS_title =  strcat('20201118 SCS Stim on/off Activation - ROI-', num2str(ns));
%         SCS_roi_hdr  = strcat('Temporal Plot:  ROI- ', num2str(ns));
%         SCS_stim_data  = eval(strcat('scs_roi_', num2str(ns)));
%         
%         scs_shft_stimData = fUS_data_shifted(SCS_stim_data, scs_stim_offset)';
%         [scs_stimON_trials_avg, scs_stimOFF_trials_avg] = get_stimONOff_trials_vec(scs_shft_stimData, stim_section_times, ONnTpts, OFFnTpts, scs_section_len);
% 
%         scs_bsline_ref = mean(mean(scs_stimOFF_trials_avg(19:29, :)));
%         scs_stimON_trials_norm = 100*(scs_stimON_trials_avg - scs_bsline_ref)/scs_bsline_ref;
%         scs_stimOFF_trials_norm = 100*(scs_stimOFF_trials_avg - scs_bsline_ref)/scs_bsline_ref;
%         
%         % *** Display Results
%         t_off_scs = 1 : size(scs_stimON_trials_avg, 1);
%         t_on_scs = size(scs_stimON_trials_avg, 1) : size(scs_stimON_trials_avg, 1) + size(scs_stimOFF_trials_avg, 1) - 1;
% 
%         figure; 
%         subplot(212); plot(SCS_stim_data); title(SCS_roi_hdr); ylabel('PD (a.u.)'); xlabel('time (s)');
%         for sm = 1 : nTrials * 2
%              if(mod(sm, 2) == 1)
%                 xoff = (sm - 1)*len_StimOFF_scs + (-1)*scs_stim_offset;
%                 if(xoff >= 0)
%                     xline(xoff,'--r',{'Off'}); 
%                 end        
%              else
%                  xon = (sm - 1)*len_StimON_scs + (-1)*scs_stim_offset;
%             if(xon >= 0)
%                 xline(xon,'--g',{'ON'});
%             end
%              end
%         end
%         
%         subplot(211); hold on;
%         shadedErrorBar(t_off_scs, scs_stimOFF_trials_norm', {@mean,@std},'lineprops','-b','transparent',1, 'patchSaturation',0.1);
%         shadedErrorBar(t_on_scs, scs_stimON_trials_norm', {@mean,@std},'lineprops','-r','transparent',1, 'patchSaturation',0.1);
%         xline(1,'--k',{'StimOFF'}); 
%         xline(len_StimON_scs,'--k',{'StimON'});
%         title(strcat('  ', SCS_title));
%         xlabel('time (s)');
%         ylabel('% Change');
%         hold off;
%         
% end

%% Auxiliary Functions 

% **** Shift Data
function fus_shift = fUS_data_shifted(rData, n_shift)
    d_len = length(rData);
  
    if(n_shift > 0)
        tpData_1 = rData(1 : d_len - n_shift);
        tpData_2 = rData((d_len - n_shift + 1) : d_len);
    
        tpShifted(1 : n_shift) = tpData_2;
        tpShifted((n_shift + 1) : (d_len)) = tpData_1;
    end
   if(n_shift < 0)
       n_shift = abs(n_shift);
        tpData_1 = rData(1 : n_shift);
        tpData_2 = rData((n_shift + 1) : d_len);
    
        tpShifted(1 : d_len - n_shift) = tpData_2;
        tpShifted((d_len - n_shift + 1) : (d_len)) = tpData_1;
   end
   if(n_shift == 0)
        tpShifted = rData';
   end
   fus_shift = tpShifted;    
end

% ****** Function returns two matrices: Given matrix of all stim ON or OFF data return trials and means of trials
% assumes equals trial epoch lengths (stimulation ON + OFF)
function [on_Stim, off_Stim] = get_stimONOff_trials_vec(stimData, stimDataTimes, on_nTpts, off_nTpts, secLen)

    stim_section_times = stimDataTimes; 
    stim_sections = 1 + floor((stim_section_times)/secLen);
    nTrials = stim_sections(end);
    
    onStim = []; offStim = [];
    for scs = 1 : nTrials
        nxt_sec = stimData(find(stim_sections == scs));

        tp_off = nxt_sec(1 : on_nTpts);
        tp_on = nxt_sec((on_nTpts + 1) : (on_nTpts + off_nTpts));

        onStim = [onStim tp_on];
        offStim = [offStim tp_off];
    end
    
    on_Stim = onStim;
    off_Stim = offStim;
end

