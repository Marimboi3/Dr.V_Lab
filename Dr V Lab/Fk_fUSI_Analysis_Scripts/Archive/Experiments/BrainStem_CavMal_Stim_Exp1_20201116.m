%% *** 20201116 GEKO Brain Stem CavMal Pre and Post Resection Stimulation Data Analysis ****

% *** Init 
clear all; close all; 
clearvars; clc;

rec_acqPeriod = 1;
rec_acqFreq = 1 / rec_acqPeriod;

% *** Load Motion Corrected Data ******
Fn_BSrec_pre_mc = '20201116_preresectionwithandwithoutGEKO_mc_1_11_0.acq';
load_BSrec_pre = load(Fn_BSrec_pre_mc, '-mat');
BSrecFull_Data_pre = squeeze(load_BSrec_pre.Acquisition.Data);   
BSrecFull_times_pre = load_BSrec_pre.Acquisition.T;

Fn_BSrec_post_mc = '20201116_postresectionwithandwithoutGEKO_mc_1_11_0.acq';
load_BSrec_post = load(Fn_BSrec_post_mc, '-mat');
BSrecFull_Data_post = squeeze(load_BSrec_post.Acquisition.Data);   
BSrecFull_times_post = load_BSrec_post.Acquisition.T;

BSrec_FullData_pre_avg = squeeze(mean(permute(BSrecFull_Data_pre, [2 1 3]), 3));
BSrec_FullData_post_avg = squeeze(mean(permute(BSrecFull_Data_post, [2 1 3]), 3));

%  *** Display Average Image ***
figure; 
subplot(121);
colormap(turbo); imagesc(BSrec_FullData_pre_avg); title('BS CavMal PreResection GekoStim fUS Average');
subplot(122); imagesc(BSrec_FullData_post_avg); title('BS CavMal PostResection GekoStim fUS Average');

% *** Load Pre/Post Resection ROIs averaged data *****
rois_BS_pre = '2020116_preresectionwithandwithoutGEKO_ROI_data_1_11_0.txt';
rois_BS_post = '2020116_postresectionwithandwithoutGEKO_ROI_data_1_11_0.txt';
load_bs_pre = importdata(rois_BS_pre);
load_bs_post = importdata(rois_BS_post);

% *** Pre Resection data *****
pre_roi_times = load_bs_pre.data(:, 1);
pre_roi_1 = load_bs_pre.data(:, 2);
pre_roi_2 = load_bs_pre.data(:, 3);
pre_roi_3 = load_bs_pre.data(:, 4);
pre_roi_4 = load_bs_pre.data(:, 5);
pre_roi_5 = load_bs_pre.data(:, 6);
pre_roi_6 = load_bs_pre.data(:, 7);
pre_roi_7 = load_bs_pre.data(:, 8);
pre_roi_8 = load_bs_pre.data(:, 9);
pre_roi_9 = load_bs_pre.data(:, 10);
pre_roi_10 = load_bs_pre.data(:, 11);
pre_roi_11 = load_bs_pre.data(:, 12);
pre_roi_12 = load_bs_pre.data(:, 13);
pre_roi_13 = load_bs_pre.data(:, 14);
pre_roi_14 = load_bs_pre.data(:, 15);
pre_roi_15 = load_bs_pre.data(:, 16);
pre_roi_16 = load_bs_pre.data(:, 17);
% pre_roi_17 = load_bs_pre.data(:, 18);
% pre_roi_18 = load_bs_pre.data(:, 19);
% pre_roi_19 = load_bs_pre.data(:, 20);
% pre_roi_20 = load_bs_pre.data(:, 21);
% pre_roi_21 = load_bs_pre.data(:, 22);
% pre_roi_22 = load_bs_pre.data(:, 23);
% pre_roi_23 = load_bs_pre.data(:, 24);
% pre_roi_24 = load_bs_pre.data(:, 25);
% pre_roi_25 = load_bs_pre.data(:, 26);
% pre_roi_26 = load_bs_pre.data(:, 27);
% pre_roi_27 = load_bs_pre.data(:, 28);

nPre_rois = 16;

% *** Post Resection data *****
post_roi_times = load_bs_post.data(:, 1);
post_roi_1 = load_bs_post.data(:, 2);
post_roi_2 = load_bs_post.data(:, 3);
post_roi_3 = load_bs_post.data(:, 4);
post_roi_4 = load_bs_post.data(:, 5);
post_roi_5 = load_bs_post.data(:, 6);
post_roi_6 = load_bs_post.data(:, 7);
post_roi_7 = load_bs_post.data(:, 8);
post_roi_8 = load_bs_post.data(:, 9);
post_roi_9 = load_bs_post.data(:, 10);
post_roi_10 = load_bs_post.data(:, 11);
post_roi_11 = load_bs_post.data(:, 12);
post_roi_12 = load_bs_post.data(:, 13);
post_roi_13 = load_bs_post.data(:, 14);
post_roi_14 = load_bs_post.data(:, 15);
post_roi_15 = load_bs_post.data(:, 16);
post_roi_16 = load_bs_post.data(:, 17);
post_roi_17 = load_bs_post.data(:, 18);
post_roi_18 = load_bs_post.data(:, 19);
post_roi_19 = load_bs_post.data(:, 20);
post_roi_20 = load_bs_post.data(:, 21);
post_roi_21 = load_bs_post.data(:, 22);
post_roi_22 = load_bs_post.data(:, 23);
post_roi_23 = load_bs_post.data(:, 24);
post_roi_24 = load_bs_post.data(:, 25);
post_roi_25 = load_bs_post.data(:, 26);
post_roi_26 = load_bs_post.data(:, 27);
post_roi_27 = load_bs_post.data(:, 28);
post_roi_28 = load_bs_post.data(:, 29);
post_roi_29 = load_bs_post.data(:, 30);
post_roi_30 = load_bs_post.data(:, 31);

nPost_rois = 30;

% *** Display  Pre Resection Results ****
   for npre = 1 : nPre_rois
       
        calMal_title_pre =  strcat(' GEKO BrainStem CavMal Stim on/off Activation (PreRes) - ROI- ', num2str(npre));
        calMal_roi_hdr_pre  = strcat('Temporal Plot:  ROI- ', num2str(npre));
        calMal_stim_data_pre  = eval(strcat('pre_roi_', num2str(npre)));
        
        % ***** CalMal Trial Stimulation ON and OFF Epoch Matices *****
        calMal_stimON_pre = [calMal_stim_data_pre(60 : 120) calMal_stim_data_pre(185 : 245)] ;
        calMal_stimOFF_pre= [calMal_stim_data_pre(127 : 181) calMal_stim_data_pre(246 : 300)];
        cmal_bsline_pre = mean(mean(calMal_stimOFF_pre((44 : 54), :)));

        calMal_stimON_norm_pre = 100*((calMal_stimON_pre - cmal_bsline_pre) / cmal_bsline_pre);
        calMal_stimOFF_norm_pre = 100*((calMal_stimOFF_pre - cmal_bsline_pre) / cmal_bsline_pre);
        cm_off_pre = 1 : 55; cm_on_pre = 55 : 115;
        
        % *** Display Results
        figure; 
        subplot(211);
        hold on;
        shadedErrorBar(cm_off_pre, calMal_stimOFF_norm_pre', {@mean,@std},'lineprops','-b','transparent',1, 'patchSaturation',0.1);
        shadedErrorBar(cm_on_pre, calMal_stimON_norm_pre', {@mean,@std},'lineprops','-r','transparent',1, 'patchSaturation',0.1);
        xline(1,'--r',{'StimOFF'}); 
        xline(54,'--g',{'StimON'});
        title(strcat('  ', calMal_title_pre));
        xlabel('time (s)');
        ylabel('% change ');
        hold off;
        subplot(212);   
        plot(calMal_stim_data_pre); title(calMal_roi_hdr_pre); ylabel('PD (a.u.)'); xlabel('time (s)');
        xline(1,'--r',{'OFF'}); xline(60,'--g',{'ON'}); xline(127,'--r',{'OFF'}); xline(185,'--g',{'ON'}); xline(245,'--r',{'OFF'});
   end
   
   % *** Display  Post Resection Results ****
%    for npos = 1 : nPost_rois
%        
%         calMal_title_post =  strcat(' GEKO BrainStem CavMal Stim on/off Activation (PostRes) - ROI- ', num2str(npos));
%         calMal_roi_hdr_post = strcat('Temporal Plot:  ROI- ', num2str(npos));
%         calMal_stim_data_post = eval(strcat('post_roi_', num2str(npos)));
%         
%         % ***** CalMal Trial Stimulation ON and OFF Epoch Matices *****
%         calMal_stimON_post = [calMal_stim_data_post(51 : 114) calMal_stim_data_post(182 : 245)] ;
%         calMal_stimOFF_post = [calMal_stim_data_post(126 : 180) calMal_stim_data_post(246 : 300)];
%         cmal_bsline_post = mean(mean(calMal_stimOFF_post((44 : 54), :)));
% 
%         calMal_stimON_norm_post = 100*((calMal_stimON_post - cmal_bsline_post) / cmal_bsline_post);
%         calMal_stimOFF_norm_post = 100*((calMal_stimOFF_post - cmal_bsline_post) / cmal_bsline_post);
%         cm_off_post = 1 : 55; cm_on_post = 55 : 118;
%         
%         % *** Display Results
%         figure; 
%         subplot(211);
%         hold on;
%         shadedErrorBar(cm_off_post, calMal_stimOFF_norm_post', {@mean,@std},'lineprops','-b','transparent',1, 'patchSaturation',0.1);
%         shadedErrorBar(cm_on_post, calMal_stimON_norm_post', {@mean,@std},'lineprops','-r','transparent',1, 'patchSaturation',0.1);
%         xline(1,'--r',{'StimOFF'}); 
%         xline(55,'--g',{'StimON'});
%         title(strcat('  ', calMal_title_post));
%         xlabel('time (s)');
%         ylabel('% change ');
%         hold off;
%         subplot(212);   
%         plot(calMal_stim_data_post); title(calMal_roi_hdr_post); ylabel('PD (a.u.)'); xlabel('time (s)');
%         xline(1,'--r',{'OFF'}); xline(51,'--g',{'ON'}); xline(125,'--r',{'OFF'}); xline(182,'--g',{'ON'}); xline(245,'--r',{'OFF'}); 
%         
%    end




