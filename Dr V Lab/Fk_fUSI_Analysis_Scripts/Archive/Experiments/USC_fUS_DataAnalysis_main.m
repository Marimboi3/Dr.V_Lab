%
% --------- USC Clinical fUS Spinal Cord/Brain data  Analysis   -----------------
%

%% Initialize and Load Files
clear all; close all; 
clearvars;
clc;

% **** Name Files to Load  *******
% FN_baseline = '20201118_baseline_1.acq';
% FN_baselineSaline = '20201118_baseline_2.acq';
% Fn_stimulation_pre = '20201116_preresectionwithandwithoutGEKO.acq';
% Fn_stimulation_post = '20201116_postresectionwithandwithoutGEKO.acq';
% Fn_stimulation_pns = '20201118_PNS Stimulation 10 minutes.acq';
% Fn_stimulation_scs = '20201118_SCS STimulation 10 minutes.acq';
% Fn_stimulation_rMVD = '20201202_RightMVD.acq';
Fn_stimulation_uds = 'UDS filling.acq';
% Fn_stimulation_bypass = 'bypass12_10.acq';


% **** Load fUS 4D Image Data Files****
Fn_1 = Fn_stimulation_uds; 
Fn_1_load = load(Fn_1, '-mat');
Fn_1_data = Fn_1_load.Acquisition.Data; 

% Fn_2 = Fn_stimulation_post; 
% Fn_2_load = load(Fn_2, '-mat');
% Fn_2_data = Fn_2_load.Acquisition.Data; 

%% Correct for Motion

 fUSdataC_1 = fUSdata_psrd_correction(Fn_1_data, 1, 9, 30, 0); % Last input: 1 = perform drift correction first, 0 = after rgid motion correction 


%% Visualize Corrected fUS Image Data

title_mc = Fn_1(1 : length(Fn_1) - 4);
vid_mc_fn = strcat(title_mc, '_MC');
imageSeqToMovie(vid_mc_fn, fUSdataC_1, size(fUSdataC_1, 4), 'turbo', strcat(title_mc, ' MC'), 0);  % 1=save movie, 0 = don't save 

%% Save Motion corrected data

MC_FileName = 'UDS_filling_mc_1_9_0.acq';
save_filtered_fUSdata(MC_FileName, fUSdataC_1, Fn_1_load);



