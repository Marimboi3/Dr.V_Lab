%%  ******** SC Epidural Electrical Stimulation and spinal cord imaging Data Analysis ****

%% Initialize workspace
clear all; close all; 
clearvars;
clc;

%% Init fUS SC data .scan file names and labeals
Fn_epidSCS_P1 = '20201118_SCS STimulation 10 minutes.acq';
Fn_epidSCS_P2  = 'SCS.acq';
Fn_epidSCS_P3  = '20210817_T10laminectomy_spinalcordstim_implant_2Dscan_2_scs_stim.scan';
Fn_epidSCS_P4  = 'scs_4.scan';
Fn_epidSCS_P5  = 'SCS_P5.scan';

FN_all = {Fn_epidSCS_P1 Fn_epidSCS_P2 Fn_epidSCS_P3 Fn_epidSCS_P4 Fn_epidSCS_P5};
FN_labels = { ' ESCS-P1'  ' ESCS-P2'  ' ESCS-P3' ' ESCS-P4' ' ESCS-P5'};

%% Filter/Pre-process all patient fUSI data 
sTart = 1;
lp1 = 0.0001;
lp2 = 0.02;
Fq = 1;
parameters = [1, 0, 0, 40, 1, [lp1, lp2, Fq], 0, 11];   % Set prprocessing parameters = [rgd p_rgd dtrn p_dtrn lwps p_lwps medn p_medn]
[epidSCS_3D_sc_raw, epidSCS_3D_sc_mc, epidSCS_times] = preProcess_epidSCData(FN_all, parameters, sTart);   %give Cells with raw and filtered 3D fUSI pD matrices

% Save data cells
% save('ESCS_3D_data_rawFiltered.mat', 'epidSCS_3D_sc_raw', 'epidSCS_3D_sc_mc', 'epidSCS_times','-append');

%% Visualize raw and filteres fUSI SC Data
Patients_contrast_raw = {[0.0 0.02] [0.0 0.05] [0.0 0.04] [0.0 0.04] [0.0 0.04]};
Patients_contrast_mc = {[0.0 0.02] [0.0 0.03] [0.0 0.04] [0.005 0.06] [0.005 0.06]};

%% Show short movie of the raw and filtered data
close all;
P_num = 4;  % patient number
mov_len = 150; % movie length (s)
next_patient = FN_labels{P_num};

% *** View raw data
title_raw = strcat(next_patient, '  Raw'); 
raw_FN = strcat(next_patient, '_raw'); 
raw_data = epidSCS_3D_sc_raw{P_num};
fUSimageToMovie_enhanced(raw_FN, raw_data, mov_len, 'hot', title_raw, 0, Patients_contrast_raw{P_num});  % 1=save movie, 0 = don't save 

% *** View filtered data
title_filt = strcat(next_patient, '  Filtered'); 
filt_FN = strcat(next_patient, '_mc'); 
filtered_data = epidSCS_3D_sc_mc{P_num};
fUSimageToMovie_enhanced(filt_FN, filtered_data, mov_len, 'hot', title_filt, 0, Patients_contrast_mc{P_num});  % 1=save movie, 0 = don't save 

%% Show 2D power Doppler mean/max  figures (Can save file)
close all;
Patients_contrast_2D = {[0.0 0.03] [0.0 0.03] [0.005 0.06] [0.008 0.1] [0.005 0.07]};

for pts = 1 : size(epidSCS_3D_sc_mc,2)
    
    mc3D = epidSCS_3D_sc_mc{pts};           %next mc data
    dataRng = 100 : size(mc3D, 3) - 100;       % data range
    
    mxAll = max(mc3D(:, :, dataRng), [], 3);    % 2D max in time dim
    mnAll = mean(mc3D, 3);                              % 2D minimum
    
    figure;
    colormap(hot);
    subplot(121); imagesc(imadjust(mnAll / max(mnAll, [], 'all'), Patients_contrast_2D{pts})); title(strcat(FN_labels{pts}, ' Mean all'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'fontsize', 14); 
    subplot(122); imagesc(imadjust(mxAll / max(mxAll, [], 'all'), Patients_contrast_2D{pts})); title(strcat(FN_labels{pts}, ' Max all'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'fontsize', 14); 
    set(gcf, 'Position',  [50, 50, 1220, 400]);
    % nextPat = strcat('epid_MxMn', FN_labels{P_num}, '.png');
    % saveas(gcf, nextPat); 
    
end

%% Data analysis
% patient protocols : Patients 1 - 3 = 10 Stim OFF - ON cycles = 600 (s) == 10 min
% Patients 4 - 5 = 3min baseline, 10 Stim ON - OFF cycles and 2min rest
acqPeriod = 1;                                                   % fUSI acquisition period (s)
acqFreq = 1 / acqPeriod;   
len_StimON_scs = 30;
len_StimOFF_scs = 30;
scs_section_len = len_StimON_scs + len_StimOFF_scs;                           % length of Stimulation ON and OFF
baseline_Ranges = {10:30 10:30 10:30 151:180 151:180 };
escsData_Ranges = {1:600 1:600 1:600 151:750 151:750 };

%% Plot global raw and filtered time series
close all;
stim_color = [1 0 1];

for pts = 1 : size(epidSCS_3D_sc_raw,2)
    
    raw3D = epidSCS_3D_sc_raw{pts};
    mc3D = epidSCS_3D_sc_mc{pts};
    next_times = epidSCS_times{pts};

    bslneRng = baseline_Ranges{pts};
    glo_raw_pct = getPercentChange_1D(squeeze(mean(mean(raw3D))), bslneRng);
    glo_mc_pct = getPercentChange_1D(squeeze(mean(mean(mc3D))), bslneRng);

    figure;
   
    set(gca,'fontsize',15); 
    title(strcat('Global Mean - ', FN_labels{pts}));
    patch([baseline_Ranges{pts}(end) baseline_Ranges{pts}(end) + 570 baseline_Ranges{pts}(end) + 570 baseline_Ranges{pts}(end)]/60, ...
        [-100 -100 100 100], stim_color, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
 hold on;
    plot(next_times/60, glo_raw_pct, '-b');
    plot(next_times/60, glo_mc_pct, '-r');

    grid;
    xline((baseline_Ranges{pts}(end) + 570)/60,'--r',{'StimOFF'}, 'fontsize',15); 
    xline(baseline_Ranges{pts}(end)/60,'--g',{'StimON'}, 'fontsize',15);
    ylabel('% pD signal change');
    xlabel('t (min)');
    yline(0,'-k',{''}); 
    ylim([-100 100]);
    legend('ESCS-ON', 'meanGlobal-Raw','meanGlobal-Filtered', 'Location','south');
    hold off;
    
end

%% Show individual trial and mean trials Stim-ON percent changes movie
close all;
p_num = 2;
Patients_contrasts = {[0.0 0.03] [0.0 0.03] [0.005 0.06] [0.008 0.1] [0.005 0.07]};
trsh_levels = {10 10 10 30 10};
cmap_Ranges = {5e1 5e1 5e1 10e1 5e1};

next_scs_3D_data = epidSCS_3D_sc_mc{p_num}(:,:, escsData_Ranges{p_num});      % load 3d trial data
scs_times = epidSCS_times{p_num}(1:600);

% ***** Get Stimulation ON and OFF Epoch Matices *****
ONnTpts = len_StimON_scs * acqFreq;
OFFnTpts = len_StimOFF_scs * acqFreq;
[scs_stimOFF_trials, scs_stimON_trials, offStmCell, onStmCell] = get_stimONOff_trials_cells(next_scs_3D_data, scs_times, ONnTpts, OFFnTpts, scs_section_len);

% *** Get mean stimulation on/off epochs ****
offStmCell_avg = (offStmCell{1} + offStmCell{2} + offStmCell{3}+ offStmCell{4} + offStmCell{5} + offStmCell{6} + offStmCell{7} + offStmCell{8} ...
    + offStmCell{9} + offStmCell{10}) / 10;
onStmCell_avg = (onStmCell{1} + onStmCell{2}+ onStmCell{3}+ onStmCell{4}+ onStmCell{5}+ onStmCell{6}+ onStmCell{7}+ onStmCell{8} ...
    + onStmCell{9} + onStmCell{10}) / 10;

% Get mean stimulation ON activation, averaged across trials
onoffStm_avg = concat3D_mat(offStmCell_avg, onStmCell_avg);     %Concat mean StimOFF and StimON
mnTBsln_rng = 24 : 30;
nextP_pchng = getPercentChange_3D(onoffStm_avg, mnTBsln_rng);

max3D_n = max(next_scs_3D_data, [], 3);
mx_mc3D_n = imadjust(max3D_n/max(max3D_n,[],'all'), Patients_contrasts{p_num});  

zSc = next_scs_3D_data(:,:,10);
wdBn_T = 2; wdBn_B = 5;
wdBn_R = 2; wdBn_L = 2;
bndH = zeros(size(zSc, 1), size(zSc, 2)); bndH(wdBn_T : size(zSc, 1) - wdBn_B, :) = 1.0;
bndV = zeros(size(zSc, 1), size(zSc, 2)); bndV(:, wdBn_L : size(zSc, 2) - wdBn_R) = 1.0;
mks_boundary = bndV.*bndH;

Guss_hft = fspecial('gaussian', 7, 0.9);
h = fspecial('disk',1.2); % disk size defined here
flt_param = [7, 1.2];
aLevel = trsh_levels{p_num};
bLevel = -trsh_levels{p_num};
cmapRng = [-cmap_Ranges{p_num} cmap_Ranges{p_num}];
colorParams = {'[0.7 0.85 0.5]' '[1.0 0.6 0.5]'};

% Show mean trial StimON Activation
nshft = -24;  % shift amount
next_section = nextP_pchng;
next_section_shift = circshift(next_section, nshft, 3);  % shift data

figure;
for im = 1 : size(next_section_shift,3)
    nxt_slc = next_section_shift(:, :, im).*mks_boundary;
    if(im < 6)
        nxtTle = strcat('\color[rgb]{0 0.5 0}Stimulation OFF ( ', num2str(im - 6), 's)');
        Stim = 0;
    end
    if(im  > 5 )
        nxtTle = strcat('\color[rgb]{0.5 0 0}Stimulation ON ( ', num2str(im - 6), 's)');

        Stim = 1;
    end
    if(im > 35)
        nxtTle = strcat('\color[rgb]{0 0.5 0}Stimulation OFF ( ', num2str(im - 6), 's)');
        Stim = 0;
    end
        
    disp(nxtTle); 
    mFrame = activation_slices(nxt_slc, mx_mc3D_n, aLevel, bLevel, flt_param, cmapRng, nxtTle, Stim, colorParams);
    pctMv(im) = mFrame;
    
end

%% Save video
%    vFn = 'Sup_Vid_1d';
%    frms = 1;
%    saveVid(vFn, pctMv, frms);

%% 
% Get single/individual Trial stimulation ON activation
nTrial = 4;
onoffStm_trial = concat3D_mat(offStmCell{nTrial+1}, onStmCell{nTrial});     %Concat mean StimOFF and StimON
TBsln_rng = 24 : 30;
nextTrial_pchng = getPercentChange_3D(onoffStm_trial, TBsln_rng);

trsh_levels = {10 25 10 30 10};
cmap_Ranges = {5e1 10e1 5e1 10e1 5e1};
aLevel = trsh_levels{p_num};
bLevel = -trsh_levels{p_num};
cmapRng = [-cmap_Ranges{p_num}  cmap_Ranges{p_num}];
colorParams = {'[0.7 0.85 0.5]' '[1.0 0.6 0.5]'};

% Show mean trial StimON Activation
nshft = -24;  % shift amount
next_section = nextTrial_pchng;
next_section_shift = circshift(next_section, nshft, 3);  % shift data

figure;
for im = 1 : size(next_section_shift,3)
    nxt_slc = next_section_shift(:, :, im).*mks_boundary;
    if(im < 6)
        nxtTle = strcat('\color[rgb]{0 0.5 0}Stimulation OFF ( ', num2str(im - 6), 's)');
        Stim = 0;
    end
    if(im  > 5 )
        nxtTle = strcat('\color[rgb]{0.5 0 0}Stimulation ON ( ', num2str(im - 6), 's)');

        Stim = 1;
    end
    if(im > 35)
        nxtTle = strcat('\color[rgb]{0 0.5 0}Stimulation OFF ( ', num2str(im - 6), 's)');
        Stim = 0;
    end
        
    disp(nxtTle); 
    mFrame = activation_slices(nxt_slc, mx_mc3D_n, aLevel, bLevel, flt_param, cmapRng, nxtTle, Stim, colorParams);
    pctMv(im) = mFrame;
    
end

%%
%    vFn = 'Sup_Vid_1d';
%    frms = 1;
%    saveVid(vFn, pctMv, frms);

%% Generate statistical parametric maps for StimON relative to "baseline" 
close all;

P_number = 4; % specify patient number
Patients_contrasts = {[0.0 0.015] [0.0 0.015] [0.009 0.025] [0.02 0.05] [0.008 0.05]};
pThresh_levels = {0.13 0.25 0.35 0.355 0.35};
nxt_Cont = Patients_contrasts{P_number};

next_scs_3D_data = epidSCS_3D_sc_mc{P_number}(:,:, escsData_Ranges{P_number});      % load 3D trial data
scs_times = epidSCS_times{P_number}(1:600);

% ***** Get Stimulation ON and OFF Epoch Matices *****
ONnTpts = len_StimON_scs * acqFreq;
OFFnTpts = len_StimOFF_scs * acqFreq;
[scs_stimOFF_trials, scs_stimON_trials, offStmCell, onStmCell] = get_stimONOff_trials_cells(next_scs_3D_data, scs_times, ONnTpts, OFFnTpts, scs_section_len);

% *** Get mean stimulation on/off epochs ****
offStmCell_avg = (offStmCell{1} + offStmCell{2} + offStmCell{3}+ offStmCell{4} + offStmCell{5} + offStmCell{6} + offStmCell{7} + offStmCell{8} ...
    + offStmCell{9} + offStmCell{10}) / 10;
onStmCell_avg = (onStmCell{1} + onStmCell{2}+ onStmCell{3}+ onStmCell{4}+ onStmCell{5}+ onStmCell{6}+ onStmCell{7}+ onStmCell{8} ...
    + onStmCell{9} + onStmCell{10}) / 10;

% Set baseline for reference
bslnRng = 24:30;
bsln_control = offStmCell_avg(:, :, bslnRng);  

% Set up SPM parameters
sigAlpha = 5e-2; % t-test significance level
fdr = 1;    % 1 means t-test with FDR correction
boundVals = [2 2 2 2]; % T, B, R, L
threshRng = [-pThresh_levels{P_number}  pThresh_levels{P_number}];
spmRng = [-3 3];  % SPM axis display range
gFilter = fspecial('gaussian', 7, 1.2);     % overlay filter
nt_title = strcat('SPM : StimON wrt StimOFF',  FN_labels{P_number}, ' (p<' , num2str(sigAlpha), ')');
spmParams = {sigAlpha fdr gFilter nt_title spmRng threshRng boundVals nxt_Cont}; % SPM papameters

figure;
next_spm = get_SPM(bsln_control, onStmCell_avg, spmParams);         % Get SPM figure
set(gcf, 'Position',  [450, 250, 800, 600]);

%% Select and plot Event related average ERA Curves from significant SPM regions
roi_select1 = roipoly(); roi_select2 = roipoly(); roi_select3 = roipoly(); roi_select4 = roipoly(); roi_select5 = roipoly(); roi_select6 = roipoly(); roi_select7 = roipoly(); 
roi_select8 = roipoly(); roi_select9 = roipoly(); roi_select10 = roipoly(); 
% roi_select11 = roipoly();  roi_select12 = roipoly(); roi_select13 = roipoly(); roi_select14 = roipoly(); 
%  roi_select15 = roipoly(); roi_select16 = roipoly();  roi_select17 = roipoly(); roi_select18 = roipoly();  roi_select19 = roipoly(); roi_select20 = roipoly(); roi_select21 = roipoly(); 
%  roi_select22 = roipoly(); roi_select23 = roipoly();roi_select24 = roipoly();
%  
mask1 = 1.0*roi_select1; mask2 = 1.0*roi_select2; mask3 = 1.0*roi_select3; mask4 = 1.0*roi_select4; mask5 = 1.0*roi_select5;  mask6 = 1.0*roi_select6; 
mask7 = 1.0*roi_select7;  mask8 = 1.0*roi_select8; mask9 = 1.0*roi_select9; mask10 = 1.0*roi_select10; 
% mask11 = 1.0*roi_select11; mask12 = 1.0*roi_select12; 
% mask13 = 1.0*roi_select13; mask14 = 1.0*roi_select14; mask15 = 1.0*roi_select15; mask16 = 1.0*roi_select16; mask17 = 1.0*roi_select17; mask18 = 1.0*roi_select18; 
% mask19 = 1.0*roi_select19; mask20 = 1.0*roi_select20; mask21 = 1.0*roi_select21; mask22 = 1.0*roi_select22; mask23 = 1.0*roi_select23; mask24 = 1.0*roi_select24;

next_scs_3D_data = epidSCS_3D_sc_mc{P_number}(:,:, escsData_Ranges{P_number});      % load 3d trial data
scs_roi_times = epidSCS_times{P_number}(1:600);
scType = FN_labels{P_number};

nOdr = 1;
scs_3D_data_ft = medfilt1(next_scs_3D_data, nOdr, [], 3); % median filter

nshft = 0;  % used to control lags
sc3D_shift = circshift(scs_3D_data_ft, nshft, 3);  % shift data
mc_3D = double(sc3D_shift);

% Get ERA curves of the ROIs, pseudo-global and global
scs_roi_1 = squeeze(mean(mean(mask1 .* mc_3D)));   %   
scs_roi_2 = squeeze(mean(mean(mask2 .* mc_3D)));   %   
scs_roi_3 = squeeze(mean(mean(mask3 .* mc_3D)));   %   
scs_roi_4 = squeeze(mean(mean(mask4 .* mc_3D)));   %  
scs_roi_5 = squeeze(mean(mean(mask5 .* mc_3D)));   %  
scs_roi_6 = squeeze(mean(mean(mask6 .* mc_3D)));   % 
scs_roi_7 = squeeze(mean(mean(mask7 .* mc_3D)));
scs_roi_8 = squeeze(mean(mean(mask8 .* mc_3D)));
scs_roi_9 = squeeze(mean(mean(mask9 .* mc_3D)));
scs_roi_10 = squeeze(mean(mean(mask10 .* mc_3D)));
% scs_roi_11 = squeeze(mean(mean(mask11 .* mc_3D)));
% scs_roi_12 = squeeze(mean(mean(mask12 .* mc_3D)));
% scs_roi_13 = squeeze(mean(mean(mask13 .* mc_3D)));
% scs_roi_14 = squeeze(mean(mean(mask14 .* mc_3D)));
% scs_roi_15 = squeeze(mean(mean(mask15 .* mc_3D)));
% scs_roi_16 = squeeze(mean(mean(mask16 .* mc_3D)));
% scs_roi_17 = squeeze(mean(mean(mask17 .* mc_3D)));
% scs_roi_18 = squeeze(mean(mean(mask18 .* mc_3D)));
% scs_roi_19 = squeeze(mean(mean(mask19 .* mc_3D)));
% scs_roi_20 = squeeze(mean(mean(mask20 .* mc_3D)));
% scs_roi_21 = squeeze(mean(mean(mask21 .* mc_3D)));
% scs_roi_22 = squeeze(mean(mean(mask22 .* mc_3D)));
% scs_roi_23 = squeeze(mean(mean(mask23 .* mc_3D)));
% scs_roi_24 = squeeze(mean(mean(mask24 .* mc_3D)));

nRois = 10;
bslinRng = 24:30;
glo_roi = 10;
scs_roi_labels = {' roi-1', ' roi-2', ' roi-3', ' roi-4', ' roi-5', ' roi-6', ' roi-7', ' roi-8', ' pSeudo',  ' Global'};
nON_pts = 30; nOFF_pts = 30; trialLen = nON_pts + nOFF_pts;
peak_times =[]; peaks = [];

for ns = 1 : nRois
    
        SCS_stim_data  = eval(strcat('scs_roi_', num2str(ns)));
        
        [stimON_trials_avg, stimOFF_trials_avg] = get_stimONOff_trials(SCS_stim_data, scs_roi_times, nON_pts, nOFF_pts, trialLen);
        stimOFFON_trials_avg = [stimOFF_trials_avg ; stimON_trials_avg];
        sc_bsline = mean(mean(stimOFFON_trials_avg(bslinRng, :))); 
        
        stimOFFON_trials_pcg = 100*(stimOFFON_trials_avg - sc_bsline)/sc_bsline;
        
        nsft = -24;  % 
        stimOFFON_trials_pcg_shift = circshift(stimOFFON_trials_pcg, nsft, 1);  % shift data
        stimOFFON_trials_pcg_mn = squeeze(mean(stimOFFON_trials_pcg_shift, 2));
        [pk, t_pk] = max(stimOFFON_trials_pcg_mn);
                        
        if(pk > 0)
            sroi = strcat('roi-', num2str(ns));
            peak_times = [peak_times t_pk];
            peaks = [peaks pk];
       
            SCS_title =  strcat('SCBV  change', ' ( ', scType, ':', scs_roi_labels{ns}, ')');
            SCS_roi_hdr  = strcat(scs_roi_labels{ns}, ' mean pD signal  time plot', ' - ', scType);
       
        % *** Display Results
            t_scs = -5:1:54;
            figure;
            set(gca,'fontsize',20); 
            hold on;
            patch_color = [1 0 1];
            patch([0 30 30 0], [-350 -350 350 350], patch_color, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
            if(ns ~= glo_roi)    
                            shadedErrorBar(t_scs, mean(stimOFFON_trials_pcg_shift',1), std(stimOFFON_trials_pcg_shift') / sqrt(size(stimOFFON_trials_pcg_shift,2)),...
                'lineprops','#0000ff','transparent',1, 'patchSaturation',0.1);
            end
            if(ns == glo_roi)            
                            shadedErrorBar(t_scs, mean(stimOFFON_trials_pcg_shift',1), std(stimOFFON_trials_pcg_shift') / sqrt(size(stimOFFON_trials_pcg_shift,2)),...
                'lineprops','#000000','transparent',1, 'patchSaturation',0.1);
            end
            xline(-5,'--g',{'StimOFF'}, 'fontsize',15); 
            xline(0,'--r',{'StimON'}, 'fontsize',15);
            xline(30,'--g',{'StimOFF'},  'fontsize',15); 
            title(strcat('  ', SCS_title));
            xline(t_pk - 6,'-.k',{strcat('pkTime :', ' ', num2str(t_pk - 5), 'sec')},  'fontsize',15); 
            yline(0,'--k',{''}); 
            yline(pk,'-.b',{strcat('Peak response :', ' ', num2str(round(pk,2)), '%')},  'fontsize',15); 
            xlabel('t (s)');
            ylabel('pD signal change (%)');
            ylim([-100 250]);
            xlim([-5 54]);
            grid;
            hold off;                     
        end
end

%% Generate PDF of Peak(Max/min) response and time peak responses distribution 
    close all;
    P_number = 4; % specify patient number
    Patients_contrasts = {[0.0 0.015] [0.0 0.015] [0.009 0.025] [0.005 0.05] [0.008 0.05]};
    pThresh_levels = {0.13 0.25 0.35 0.355 0.35};
    nxt_Cont = Patients_contrasts{P_number};
    alv = pThresh_levels{P_number}; 
    blv = -pThresh_levels{P_number}; 
    colormap_Rng = {50 50 50 150};
    mxThresh = 400;
    scType = FN_labels{P_number};
    mxCaxis = [-colormap_Rng{P_number} colormap_Rng{P_number}];

    next_scs_3D_data = epidSCS_3D_sc_mc{P_number}(:,:, escsData_Ranges{P_number});      % load 3D trial data
    scs_times = epidSCS_times{P_number}(1:600);

    mxRng = 10 : 590;
    mxAll_bg = max(next_scs_3D_data(:, :, mxRng), [], 3);
    mxAll_bg_n = imadjust((mxAll_bg / max(mxAll_bg, [], 'all')), nxt_Cont);

    % ***** Get Stimulation ON and OFF Epoch Matices *****
    ONnTpts = len_StimON_scs * acqFreq;
    OFFnTpts = len_StimOFF_scs * acqFreq;
    [scs_stimOFF_trials, scs_stimON_trials, offStmCell, onStmCell] = get_stimONOff_trials_cells(next_scs_3D_data, scs_times, ONnTpts, OFFnTpts, scs_section_len);

    % *** Get mean stimulation on/off epochs ****
    offStmCell_avg = (offStmCell{1} + offStmCell{2} + offStmCell{3}+ offStmCell{4} + offStmCell{5} + offStmCell{6} + offStmCell{7} + offStmCell{8} ...
        + offStmCell{9} + offStmCell{10}) / 10;
    onStmCell_avg = (onStmCell{1} + onStmCell{2}+ onStmCell{3}+ onStmCell{4}+ onStmCell{5}+ onStmCell{6}+ onStmCell{7}+ onStmCell{8} ...
        + onStmCell{9} + onStmCell{10}) / 10;

    scs_OffOn_avg = concat3D_mat(onStmCell_avg, offStmCell_avg);
    nshft = 6;  % shift amount
    scs_OffOn_avg_shift = circshift(scs_OffOn_avg, nshft, 3);  % shift data

    bslnRng = 24:30;
    bsln_control = offStmCell_avg(:, :, bslnRng); bsln_control_mn = mean(bsln_control, 3);

    nBsln = 1 : 5; nDataRng = 6 : 60;
    nwData_3D = scs_OffOn_avg_shift(:, :, nDataRng);
    avgPct_Chg3D = double(100*((nwData_3D - bsln_control_mn) ./ bsln_control_mn));

    sigAlpha = 5e-2;
    fdr = 1;
    [pTrsh, trshIndx, zSc] = get_sigIndex_fdr(bsln_control, onStmCell_avg, sigAlpha, fdr);

    tmRng = [-5 55];

    wdBn_T = 2; wdBn_B = 2;
    wdBn_R = 2; wdBn_L = 2;
    bndH = zeros(size(zSc, 1), size(zSc, 2)); bndH(wdBn_T : size(zSc, 1) - wdBn_B, :) = 1.0;
    bndV = zeros(size(zSc, 1), size(zSc, 2)); bndV(:, wdBn_L : size(zSc, 2) - wdBn_R) = 1.0;
    mks_boundary = bndV.*bndH;

    % zSc_trsh = double(zSc.*trshIndx);
    zSc_trsh = double(zSc.*trshIndx.*mks_boundary);   %Without the border
    zSc_trsh_abvIndx = (zSc_trsh > alv)*1.0;  
    zSc_trsh_blwIndx = (zSc_trsh < blv)*1.0; 

    zSc_abv_gus = imfilter(zSc.*zSc_trsh_abvIndx, Guss_hft, 'replicate');
    zSc_blw_gus = imfilter(zSc.*zSc_trsh_blwIndx, Guss_hft, 'replicate');


    [mxVal, t_mxVal] = max(avgPct_Chg3D, [], 3);
    [mnVal, t_mnVal] = min(avgPct_Chg3D, [], 3);

    mxVal_gus = imfilter(mxVal, Guss_hft, 'replicate');    %mxVal_gus = mxVal_gus.*mks_boundary;
    mxVal_sig = mxVal_gus.*((zSc_abv_gus > alv)*1.0);

    mnVal_gus = imfilter(mnVal, Guss_hft, 'replicate');     %mnVal_gus = mnVal_gus.*mks_boundary;
    mnVal_sig = mnVal_gus.*((zSc_blw_gus < blv)*1.0);

    t_Guss_hft = fspecial('gaussian', 7, 0.9);
    t_mxVal_gus = imfilter(t_mxVal, t_Guss_hft, 'replicate'); 
    t_mxVal_sig = t_mxVal_gus.*((zSc_abv_gus > alv)*1.0);

    t_mnVal_gus = imfilter(t_mnVal, t_Guss_hft, 'replicate');
    t_mnVal_sig = t_mnVal_gus.*((zSc_blw_gus < blv)*1.0);

    mxVal_sigAct = mxVal.*zSc_trsh_abvIndx;
    mxVal_vec = reshape(mxVal_sigAct, size(mxVal_sigAct,1)*size(mxVal_sigAct,2), 1);
    mxVal_sub_mz = mxVal_vec(mxVal_vec ~= 0);
    mxVal_sub = mxVal_sub_mz(mxVal_sub_mz < mxThresh);
    min_mxValSub = min(mxVal_sub);

    mnVal_sigAct = mnVal.*zSc_trsh_blwIndx;
    mnVal_vec = reshape(mnVal_sigAct, size(mnVal_sigAct,1)*size(mnVal_sigAct,2), 1);
    mnVal_sub = mnVal_vec(mnVal_vec ~= 0);
    max_mnValSub = max(mnVal_sub);

    MnMxVal_vec = [mnVal_sub; mxVal_sub];

    t_mxVal_sigAct = t_mxVal.*zSc_trsh_abvIndx;
    t_mxVal_vec = reshape(t_mxVal_sigAct, size(t_mxVal_sigAct,1)*size(t_mxVal_sigAct,2), 1);  
    t_mxVal_vec = t_mxVal_vec(t_mxVal_vec ~= 0);

    t_mnVal_sigAct = t_mnVal.*zSc_trsh_blwIndx;
    t_mnVal_vec = reshape(t_mnVal_sigAct, size(t_mnVal_sigAct,1)*size(t_mnVal_sigAct,2), 1);  
    t_mnVal_vec = t_mnVal_vec(t_mnVal_vec ~= 0);

    t_MnMxVal_sub = [t_mnVal_vec; t_mxVal_vec];

% max min values overlay
    figure;
    ax1 = axes;
    imagesc(ax1, mxAll_bg_n); set(gca,'fontsize',15); 
    title(strcat('Peak pD response-changes (', scType, ')'));
    colormap(ax1,'gray');
    set(gca,'xtick',[]); set(gca,'xticklabel',[]);
    set(gca,'ytick',[]); set(gca,'yticklabel',[]);

    ax2 = axes;
    ia = imagesc(ax2, mxVal_sig, 'alphadata', mxVal_sig > 0.0);
    ct = colorbar('FontSize', 20,'Location', 'westoutside');
    ct.Label.String = 'Peak response (%)';
    colormap(ax2, jet);         
    ax2.Visible = 'off';
    caxis(ax2, mxCaxis);

    ax3 = axes;
    ib = imagesc(ax3, mnVal_sig, 'alphadata', mnVal_sig < 0.0);
% ct = colorbar('FontSize', 20,'Location', 'eastoutside');
    colormap(ax3,  jet);         
    ax3.Visible = 'off';
    caxis(ax3, mxCaxis);
    
    linkprop([ax1, ax2, ax3],'Position'); 
    MinMaxVals = getframe(gcf);
    set(gcf, 'Position',  [450, 250, 800, 600]);

% times to max overlay
    figure;
    ax1 = axes;
    imagesc(ax1, mxAll_bg_n); set(gca,'fontsize',15); 
    title(strcat('Peak reponse times (', scType, ')'));
    colormap(ax1,'gray');
    set(gca,'xtick',[]); set(gca,'xticklabel',[]);
    set(gca,'ytick',[]); set(gca,'yticklabel',[]);

    ax2 = axes;
    imagesc(ax2, t_mxVal_sig, 'alphadata', t_mxVal_sig > 0.0);
    ct = colorbar('FontSize', 20,'Location', 'westoutside');
    ct.Label.String = 'Peak time (s)';
    colormap(ax2, jet);         
    ax2.Visible = 'off';
    caxis(ax2, tmRng) ;

    ax3 = axes;
    imagesc(ax3, t_mnVal_sig, 'alphadata', t_mnVal_sig > 0.0);
    colormap(ax3,  jet);         
    ax3.Visible = 'off';
    caxis(ax3, tmRng);

    linkprop([ax1, ax2, ax3],'Position'); 
    tMinMax = getframe(gcf);
    set(gcf, 'Position',  [450, 250, 800, 600]);

    patch_color = [0 0 0];           
    % plot PDF histograms combined
     figure;
    set(gca,'fontsize',15); 
    subplot(121); set(gca,'fontsize',15); 
    hold on;
    ht1 = histogram(t_mxVal_vec, 'Normalization','Probability', 'EdgeColor', [0.5 0 0]);
    ht1.FaceColor = [1 0.1 0];
    alpha(ht1, 0.2);
    
    ht2 = histogram(t_mnVal_vec, 'Normalization','Probability', 'EdgeColor', [0 0 0.5]);
    ht2.FaceColor = [0 0.1 1];
    alpha(ht2, 0.2);
    
    mxHis = 1.05*(max([max(ht1.Values) max(ht2.Values)]));
    patch([0 30 30 0], [-1 -1 1 1], patch_color, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    xline(0,'--r',{'StimON'}, 'fontsize',12);
    xline(30,'--k',{'StimOFF'}, 'fontsize',12);
    title(strcat('PDF Peak times (', scType, ')'));
    xlabel('Peak times (s)');
    ylabel('Probability');
    xlim(tmRng);
    ylim([0 mxHis]);
    grid;
    legend('TMax empirical PDF', 'TMin empirical PDF');
    hold off;
    
    subplot(122); set(gca,'fontsize',15); 
    hold on;
    hp1 = histogram(mxVal_sub, 'Normalization','Probability');
    hp1.FaceColor = [1 0 0];
    
    hp2 = histogram(mnVal_sub, 'Normalization','Probability');
    hp2.FaceColor = [0 0 1];
    
    title(strcat('PDF Peak responses (', scType, ')'));
    grid;
    ylabel('Probability');
    xlabel('Peak signal change (%)');
    legend('Peak (max) empirical PDF', 'Peak (min) empirical PDF');
    hold off;
    set(gcf, 'Position',  [450, 250, 1200, 500]);
    
 %% Auxiliary functions
 
 function [raw3D, mc3D, scsTimes] = preProcess_epidSCData(dataCell, params, strt)
    mc_3D_files = {}; raw_3D_files = {}; scs_times = {};
       
    for a = strt : size(dataCell, 2)
        nextFN = dataCell{a};
        if(a <= 2)
            Fn_load = load(nextFN, '-mat');
            nextData = Fn_load.Acquisition.Data;
            next_times = Fn_load.Acquisition.T;
        end
        if(a > 2)
            nextData = h5read(nextFN, '/Data'); 
            next_times = h5read(nextFN, '/acqMetaData/time'); 
        end        
        nextData_mc = fUSdata_psrd_filter(nextData, params, 0); % Last input: 1 = perform drift correction first, 0 = after rgid motion correction 
        nextData_raw_perm = permute(squeeze(nextData), [2 1 3]);
        nextData_mc_perm = permute(squeeze(nextData_mc), [2 1 3]);
        
        raw_3D_files{a} = double(nextData_raw_perm);
        mc_3D_files{a} = double(nextData_mc_perm);
        scs_times{a} = double(next_times);
        
    end
    
    raw3D = raw_3D_files;
    mc3D = mc_3D_files;
    scsTimes = scs_times;
    
 end
 
 function pct = getPercentChange_1D(dat, brng)

    bsln = mean(dat(brng));
    ypct = 100*((dat - bsln) / bsln);
    
    pct = ypct;
end
 
 function mdFilt = median_filter_2D(dat, odr)
    tpNdf = [];
    for i = 1 : size(dat, 2)
           tpNdf(:, i) = medfilt1(dat(:, i), odr);
    end
    mdFilt = tpNdf;
end

function pct = getPercentChange_2D(dat, brng)

    bsln = mean(dat(brng, :), 1);
    ypct = 100*((dat - bsln) ./ bsln);    
    pct = ypct;
    
end

function mdFilt = median_filter_1D(dat, odr)
    ytp = medfilt1(dat, odr);
    mdFilt = ytp;
end

function pct = getPercentChange_3D(dat, brng)

    bsln = mean(dat(:, :, brng), 3);
    ypct = 100*((dat - bsln) ./ bsln);
    
    pct = ypct;
    
end

function frame = activation_slices(dat, bgimg, aLevel, bLevel, flt_param, cmapRng, Tle, Stim, cPm)
        
         Guss_hft = fspecial('gaussian', flt_param(1), flt_param(2));
         K_gus = imfilter(dat, Guss_hft, 'replicate');
         clf;     
         ax1 = axes;
         imagesc(ax1, bgimg);
         colormap(ax1,'gray');
         set(gca,'fontsize',15);
         set(gca,'xtick',[]); set(gca,'xticklabel',[]);
         set(gca,'ytick',[]); set(gca,'yticklabel',[]);
         xlabel( '|< --------------------------      width (12.8 mm)      -------------------------- >|');
         ylabel( '|< ---------------     depth (10.0 mm)     --------------- >|');

        if(Stim == 0)
            title(Tle,'fontsize', 20, 'BackgroundColor', cPm{1});
        end
        if(Stim == 1)
            title(Tle, 'fontsize', 20, 'BackgroundColor', cPm{2});
        end
         
         ax2 = axes;
         ima = imagesc(ax2, K_gus, 'AlphaData', K_gus > aLevel);
         
         colormap(ax2,  jet);         
         ax2.Visible = 'off';
         caxis(ax2, cmapRng);
          set(gca,'fontsize',18);

         ax3 = axes;
         imb = imagesc(ax3, K_gus, 'alphadata', K_gus < bLevel);
         colormap(ax3,  jet);         
         ax3.Visible = 'off';
         caxis(ax3, cmapRng);
          set(gca,'xtick',[]); set(gca,'xticklabel',[]);
         set(gca,'ytick',[]); set(gca,'yticklabel',[]);
          set(gca,'fontsize',15); 
         
         linkprop([ax1, ax2, ax3],'Position'); 
         set(gcf, 'Position',  [50, 50, 850, 600]);
         frame = getframe(gcf);

end
   
   
% ****** Function returns two matrices: Given matrix of all stim ON or OFF data return trials and means of trials
% assumes equals trial epoch lengths (stimulation ON + OFF)
function [off_Stim, on_Stim, offC, onC] = get_stimONOff_trials_cells(stimData, stimDataTimes, on_nTpts, off_nTpts, secLen)

    stim_section_times = stimDataTimes; 
    stim_sections = 1 + floor((stim_section_times)/secLen);
    nTrials = stim_sections(end);
    
    onStim = []; offStim = [];
    for scs = 1 : nTrials
        nxt_sec = stimData(:, :,  find(stim_sections == scs));

        tp_off = nxt_sec(:,:, 1 : on_nTpts);
        tp_on = nxt_sec(:, :, (on_nTpts + 1) : (on_nTpts + off_nTpts));

        offStim = [offStim, permute(tp_off, [1 3 2])];
        onStim = [onStim, permute(tp_on, [1 3 2])];
        onCell{scs} = tp_on;
        offCell{scs} = tp_off;
    end
    onC = onCell;
    offC = offCell;
    off_Stim = permute(offStim, [1 3 2]);
    on_Stim = permute(onStim, [1 3 2]);
    
end

function combo3D = concat3D_mat(mat1, mat2)

    tp = [permute(mat1, [1 3 2]), permute(mat2, [1 3 2])];
    combo3D = permute(tp, [1 3 2]);
    
end

   
function saveVid(fname, mov, nfrm)

        vidTp = VideoWriter(fname);
        vidTp.FrameRate = nfrm;
        open(vidTp);
        writeVideo(vidTp, mov);
        close(vidTp);
end

function spmImg = get_SPM(bslineDat, relDat, params)

    data_max = max(bslineDat, [], 3);
    sc_BgImg = imadjust(data_max / max(data_max, [], 'all'), params{8});

    [pTrsh, trshIndx, zSc] = get_sigIndex_fdr(bslineDat, relDat, params{1}, params{2});
    
    wdBn_T = params{7}(1); wdBn_B = params{7}(2);
    wdBn_R = params{7}(3); wdBn_L = params{7}(4);
    bndH = zeros(size(zSc, 1), size(zSc, 2)); bndH(wdBn_T : size(zSc, 1) - wdBn_B, :) = 1.0;
    bndV = zeros(size(zSc, 1), size(zSc, 2)); bndV(:, wdBn_L : size(zSc, 2) - wdBn_R) = 1.0;
    mks_boundary = bndV.*bndH;

    zSc_trsh = double(zSc.*trshIndx.*mks_boundary);   %With the border
    zSc_trsh_abvIndx = (zSc_trsh > params{6}(2))*1.0;  
    zSc_trsh_blwIndx = (zSc_trsh < params{6}(1))*1.0; 

    zSc_abv_gus = imfilter(zSc.*zSc_trsh_abvIndx, params{3}, 'replicate');
    zSc_blw_gus = imfilter(zSc.*zSc_trsh_blwIndx, params{3}, 'replicate');
    
    % get SPM figure
    ax1 = axes;
    imagesc(ax1, sc_BgImg); set(gca,'fontsize',15); 
    title(params{4});
    colormap(ax1,'gray');
    set(gca,'xtick',[]); set(gca,'xticklabel',[]);
    set(gca,'ytick',[]); set(gca,'yticklabel',[]);

    ax2 = axes;
    ima = imagesc(ax2, zSc_abv_gus, 'alphadata', zSc_abv_gus > params{6}(2));
    colormap(ax2,  jet);         
    ax2.Visible = 'off';
    caxis(ax2, params{5});
    c2 = colorbar('FontSize', 20,'Location', 'westoutside');
    c2.Label.String = 'z-scores';

    ax3 = axes;
    imb = imagesc(ax3, zSc_blw_gus, 'alphadata', zSc_blw_gus < params{6}(1));
    colormap(ax3,  jet);         
    ax3.Visible = 'off';
    caxis(ax3, params{5});

    linkprop([ax1, ax2, ax3],'Position');   
    spmImg = getframe(gcf);
end

function [pVals, treshIndx, zSc] = get_sigIndex_fdr(ref_mat, mat2, alpha, fdr)

            mn_mt1 = mean(ref_mat,3); 
            mn_mt2 = mean(mat2, 3);

            [rject, p_vals] = ttest2(permute(ref_mat, [3 2 1]), permute(mat2, [3 2 1]));
            p_vals = squeeze(p_vals)';

            if(fdr == 1)
                    pval_vec = reshape(p_vals, size(p_vals, 1)*size(p_vals, 2), 1);
                    [c_pvalues, c_alpha, h, extra] = fdr_BH(pval_vec, alpha);
                    pval_fdr_add = c_pvalues;
                    pval_fdr = reshape(pval_fdr_add, size(p_vals, 1), size(p_vals, 2));
                    p_vals = pval_fdr;
            end    
            sig_indx = p_vals < alpha;

            ztemp = zeros(size(mat2,1), size(mat2,2)); 
            p_tresh = ztemp;  

            p_tresh(sig_indx) = p_vals(sig_indx);
            mn_tresh = mn_mt2 - mn_mt1;

            tp_vec = reshape(mn_tresh, size(mn_tresh,1)*size(mn_tresh,2), 1);
            zSc_vec = zscore(tp_vec);
            sigZscore = reshape(zSc_vec, size(mn_tresh,1), size(mn_tresh,2));

            treshIndx = sig_indx;
            pVals = p_tresh;
            zSc = sigZscore;
end

% ****** Function returns two matrices: Given matrix of all stim ON or OFF data return trials and means of trials
% assumes equals trial epoch lengths (stimulation ON + OFF)
function [on_Stim, off_Stim] = get_stimONOff_trials(stimData, stimDataTimes, on_nTpts, off_nTpts, secLen)

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
