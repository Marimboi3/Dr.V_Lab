% --------- Sample fUSI data pre-processing pipeline - animal or Human data -----------------

%% Initialize workspace
clear; close all; clearvars; clc;

%% Init subject fUSI data .scan file names and labels

% load sample data
FN_Sub1_F1 = '20220517_transcu_SC_P1.scan'; 
FN_Sub1_F2 = 'scs_uds_82322_transcutaneous_P2.scan'; 
Sub1_FNs = {FN_Sub1_F1 FN_Sub1_F2};

FN_Sub2_F1 = 'mk_mse4_c1_5.scan'; 
FN_Sub2_F2 = 'mk_mse5_c1_5.scan'; 
Sub2_FNs = {FN_Sub2_F1 FN_Sub2_F2};

% All subjects file names and labels
all_subjects_FNs = {Sub1_FNs Sub2_FNs};
all_subjects_FLabels = {{ ' sub1-F1'  ' sub1-F2'} {' sub2-F1'  ' sub2-F2'}};

%% Filter fUSI data  
sTart = 1;
lp1 = 0.0001;
lp2 = 0.04;
Fq = 1;

refLen = 40; binwdt = 120; 
parameters = {{0 refLen binwdt} {1 lp1 lp2 Fq}};   % Set preprocessing parameters = cell1 - NormCorre cell2 - lowpass filter
exp_subjectData_all = preProcess_AnimalHumanData(all_subjects_FNs, parameters, sTart);   % give Cells with each patients ...
% exp_subjectData_all = preProcess_AnimalHumanData(all_patients_FNs, parameters, sTart, bsLnRng);   % use if normCorre image template is spgive Cells with each patients ...

% Save data cells
% save('subject_rawAndFiltered_data.mat', 'exp_subjectData_all', '-append');

%% Data analysis
acqPeriod = 1;                                                   % fUSI acquisition period (s)
acqFreq = 1 / acqPeriod;   

%% Visualize raw and corrected fUS Data
% Image parameters
subj_contrast_raw = {{[0.01 0.1] [0.0 0.015]} {[0.0 0.12] [0.0 0.2]}};
subj_contrast_mc = {{[0.01 0.2] [0.0 0.015]} {[0.01 0.2] [0.02 0.2]}};
subj_hisAdj_Levels = {{0.01 0.02} {0.01 0.02}};
subJ_bsline_ranges = {{1:120 1:120} {241:300 241:300}};

%% Visualize raw and corrected fUS Data
close all;
sub_num = 2;
FL_num = 2;

subj_3D_raw = exp_subjectData_all{sub_num}{1};
subj_3D_mc = exp_subjectData_all{sub_num}{2};
subj_3D_times = exp_subjectData_all{sub_num}{3};

nextRaw_cell = subj_3D_raw{FL_num};
nextMc_cell = subj_3D_mc{FL_num};

raw_FN = strcat(all_subjects_FLabels{sub_num}{FL_num}, '_raw');
filt_FN= strcat(all_subjects_FLabels{sub_num}{FL_num}, '_mc');
FTitles = { strcat('Vascular map pD movie  (', all_subjects_FLabels{sub_num}{FL_num}, '-raw)' )  strcat('vascular map pD movie  (', all_subjects_FLabels{sub_num}{FL_num}, '-filtered)' )};
mov_len = 100; 

% *** View raw data
title_raw = FTitles{1}; 
raw_data = nextRaw_cell;
fUSimageToMovie_enhanced(raw_FN, raw_data, mov_len, 'hot', title_raw, 0, subj_contrast_raw{sub_num}{FL_num});  % 1=save movie, 0 = don't save 

% *** View filtered data
title_filt = FTitles{2}; 
filtered_data = nextMc_cell;
fUSimageToMovie_enhanced(filt_FN, filtered_data, mov_len, 'hot', title_filt, 0,  subj_contrast_mc{sub_num}{FL_num});  % 1=save movie, 0 = don't save

%% 2D mean/max power doppler figures
close all

sub_num = 1;
mnMxStd_contrast = {{[0.025 0.8] [0.035 0.9] [0.0 0.5]} {[0.0 0.45] [0.0 0.4] [0.0 0.4]}};

nxt_subj_data_mc = exp_subjectData_all{sub_num}{2};
nxt_subj_labels = all_subjects_FLabels{sub_num};

for Fnum = 1 : size(nxt_subj_data_mc, 2)

    nxtRec = nxt_subj_data_mc{Fnum}(:, :, subJ_bsline_ranges{Fnum}{1});
    nxtLabel = nxt_subj_labels{Fnum}; 
    nxtHistLvl = subj_hisAdj_Levels{sub_num}{Fnum};

    nxtMx = max(nxtRec, [], 3);      nxtMx_n = nxtMx / max(nxtMx, [], 'all');
    nxtMn = mean(nxtRec, 3);         nxtMn_n = nxtMn / max(nxtMn, [], 'all');  
    nxtStd = std(nxtRec, 0, 3);        nxtStd_n = nxtStd / max(nxtStd, [], 'all');

    nxMx_cont = imadjust(nxtMx_n, mnMxStd_contrast{Fnum}{1}, []);
    nxMn_cont = imadjust(nxtMn_n, mnMxStd_contrast{Fnum}{2}, []);
    nxStd_cont = imadjust(nxtStd_n, mnMxStd_contrast{Fnum}{3}, []);

    nxMx_adj = adapthisteq(nxMx_cont,'clipLimit', nxtHistLvl,'Distribution','rayleigh');
    nxMn_adj = adapthisteq(nxMn_cont,'clipLimit', nxtHistLvl,'Distribution','rayleigh');
    nxStd_adj = adapthisteq(nxStd_cont,'clipLimit', nxtHistLvl,'Distribution','rayleigh');  

    figure;
    colormap(hot);
    subplot(1, 3, 1); 
    imagesc(nxMx_adj);  title(strcat(nxtLabel, ' - pD Max'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'fontsize', 9); 

    subplot(1, 3, 2); 
    imagesc(nxMn_adj);  title(strcat(nxtLabel, ' - pD Mean'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'fontsize', 9); 

    subplot(1, 3, 3); 
    imagesc(nxStd_adj);  title(strcat(nxtLabel, ' - pD Std'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'fontsize', 9); 


    set(gcf, 'Position',  [50, 150, 1350, 300]);
    nextGrp = strcat('MxMnStd_VascularMap_trancuStim_' , '.png');
%     saveas(gcf, nextGrp);

end

%% Plot global raw and filtered time series
close all;
sub_num = 2;

nxtSubj_data_raw = exp_subjectData_all{sub_num}{1};
nxtSubj_data_mc = exp_subjectData_all{sub_num}{2};
nxtSubj_times = exp_subjectData_all{sub_num}{3};
nxtSubj_labels = all_subjects_FLabels{sub_num};
nxtBsln_Rng = subJ_bsline_ranges{sub_num};

startFile = 1;

bsln_color = [0 0 0.1];
stim_color = [1 0 1];

for Fns = startFile : size(nxtSubj_data_mc,2)
    
    raw3D = nxtSubj_data_raw{Fns};
    mc3D = nxtSubj_data_mc{Fns};
    next_times = nxtSubj_times{Fns};

    bslneRng = nxtBsln_Rng{Fns};
    glo_raw_pct = getPercentChange_1D(squeeze(mean(mean(raw3D))), bslneRng);
    glo_mc_pct = getPercentChange_1D(squeeze(mean(mean(mc3D))), bslneRng);

    figure;   
    set(gca,'fontsize',20); 
    title(strcat('Global Mean - ', nxtSubj_labels{Fns}));           
    
%         patch([patient_protocols{pts}{1}(end) patient_protocols{pts}{6}(end) patient_protocols{pts}{6}(end) patient_protocols{pts}{1}(end)]/60, ...
%         [-200 -200 200 200], stim_color, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    
 hold on;
    plot(next_times/60, glo_raw_pct, '-b');
    plot(next_times/60, glo_mc_pct, '-r', LineWidth=2.0);
    patch([bslneRng(1) bslneRng(end) bslneRng(end) bslneRng(1)]/60, ...
        [-200 -200 200 200], bsln_color, 'FaceAlpha', 0.05, 'EdgeColor', 'none');

    xlb = xline(bslneRng(1)/60,'--k',{'baseline'}, 'fontsize',15); 
    xlb.LabelHorizontalAlignment = 'center';

    xls = xline(bslneRng(end)/60,'--r',{'Start-task'}, 'fontsize',15);
    xls.LabelVerticalAlignment = 'bottom';
    ylabel('% pD signal change');
    xlabel('t (min)');
    yline(0,'-k',{''}); 
    ylim([-60 60]);
    legend( 'meanGlobal-Raw','meanGlobal-Filtered', 'Location','best');
    hold off;
    
    set(gcf, 'Position',  [50, 50, 900, 550]);    
end

%% Visualize spectograph of data across pixels over time
close all;
sub_num = 1;
F_num = 1;

nxtSubj_data_mc = exp_subjectData_all{sub_num}{2};
nxtSubj_labels = all_subjects_FLabels{sub_num};
nxtBsln_Rng = subJ_bsline_ranges{sub_num};

subj_contrast_mc = {{[0.0 0.5] [0.1 0.2]} {[0.01 0.2] [0.0 0.25]}};

 mc3D_sigTrial = nxtSubj_data_mc{F_num};
 sigTrial_bslneRng = subJ_bsline_ranges{sub_num}{F_num};
 sigTrial_mc_3D_pct = getPercentChange_3D(mc3D_sigTrial, sigTrial_bslneRng);

 sigTrial_mc_3D_pct_vec = reshape(sigTrial_mc_3D_pct, size(sigTrial_mc_3D_pct,1)*size(sigTrial_mc_3D_pct,2), size(sigTrial_mc_3D_pct,3));

sigTrial_mc_3D_pct_vec_n = sigTrial_mc_3D_pct_vec/max(sigTrial_mc_3D_pct_vec,[],'all');  
sigTrial_mc_3D_pct_vec_n_adj = adapthisteq(sigTrial_mc_3D_pct_vec_n,'clipLimit', 0.025,'Distribution','rayleigh');    
sigTrial_mc_3D_pct_vec_n_adj = imadjust(sigTrial_mc_3D_pct_vec_n_adj, subj_contrast_mc{sub_num}{F_num}, []);

 figure;
 colormap(jet);
 imagesc(sigTrial_mc_3D_pct_vec_n_adj);
 title(strcat('Pixel pD change spectrograph (P-', num2str(F_num), ')'));
 xlabel('t (s)');
 ylabel('pixels');
 set(gca,'fontsize', 15); 
 set(gcf, 'Position',  [50, 50, 1300, 730]);  

%% auxiliary function

 function pct = getPercentChange_1D(dat, brng)

    bsln = mean(dat(brng));
    ypct = 100*((dat - bsln) / bsln);
    
    pct = ypct;
end


function pct = getPercentChange_3D(dat, brng)

    bsln = mean(dat(:, :, brng), 3);
    ypct = 100*((dat - bsln) ./ bsln);
    
    pct = ypct;
    
end