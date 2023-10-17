
% --------- Pre-clinical fUSI Mice MK801 dose dependant data pre-processing analysis   -----------------

%% Clear workspaces  
clear; close all; clearvars; clc;

%% Assign variables to mice pD scan files

Fn_mse_1 = 'mk_mse2_c0_1.scan';
Fn_mse_5 = 'mk_mse5_c1_5.scan';

mice_FN_all = {{Fn_mse_1 Fn_mse_5}};

baseline_Ranges_all = {241:300 241:300};

%% Load DBS data files Data anaysis fitered data


%% Pre-process all patient fUSI data 
sTart = 1;           % Set filter parameters
lp1 = 0.0001;
lp2 = 0.03;
Fq = 1;

refLen = 30; binwdt = 120; bsLnRng = 113:119;
parameters = {{1 refLen binwdt} {1 lp1 lp2 Fq}};   % Set preprocessing parameters = cell1 - NormCorre cell2 - lowpass filter
miceData_all = preProcess_mkDoseData(mice_FN_all, parameters, sTart);   % give Cells with each patients ...
% miceData_all = preProcess_mkDoseData(mice_FN_all, parameters, sTart, bsLnRng);   % use if normCorre image template is spgive Cells with each patients ...

%% Visualize raw and filteres fUSI SC Data

mse_contrast_raw = {[0.0 0.4] [0.0 0.2] };
mse_contrast_mc = {[0.0 0.5] [0.0 0.4]  };

mice_contrast_raw = {mse_contrast_raw};
mice_contrast_mc = {mse_contrast_mc};

%% Show short movie of the raw and filtered data
close all;

Exp_labels = {{'mse-2' 'mse-5'}};
m_num = 1;
nextMse_data = miceData_all{m_num};

File_num = 2;  % patient number
mov_len = 200; % movie length (s)
next_mse = Exp_labels{m_num}{File_num};

% *** View raw data
title_raw = strcat(next_mse, '  Raw'); 
raw_FN = strcat(next_mse, '_raw'); 
raw_data = nextMse_data{1}{File_num};
fUSimageToMovie_enhanced(raw_FN, raw_data, mov_len, 'hot', title_raw, 0, mice_contrast_raw{m_num}{File_num});  % 1=save movie, 0 = don't save 

% *** View filtered data
title_filt = strcat(next_mse, '  Filtered'); 
filt_FN = strcat(next_mse, '_mcFrq'); 
filtered_data = nextMse_data{2}{File_num};
fUSimageToMovie_enhanced(filt_FN, filtered_data, mov_len, 'hot', title_filt, 0, mice_contrast_mc{m_num}{File_num});  % 1=save movie, 0 = don't save 

%% Show 2D power Doppler mean/max  figures (Can save file)
close all;
m_num = 1;
nextMse_data_mc = miceData_all{m_num}{2};
mice_contrast_2D = {[0.0 0.45] [0.0 0.5] };

for pts = 1 : size(nextMse_data_mc,2)
    
    mc3D = nextMse_data_mc{pts};           %next mc data
    dataRng = baseline_Ranges_all{m_num}; %25 : size(mc3D, 3) - 25;       %  bdata range
    
    mxAll = max(mc3D(:, :, dataRng), [], 3);    % 2D max in time dim
    mnAll = mean(mc3D, 3);                              % 2D minimum
    stdAll = std(mc3D, 0, 3);                              % 2D minimum
    
    figure;
    colormap(hot);
    subplot(131); imagesc(imadjust(mnAll / max(mnAll, [], 'all'), mice_contrast_2D{pts})); title(strcat(Exp_labels{m_num}{pts}, ' bsline-Mean'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'fontsize', 14); 
    
    subplot(132); imagesc(imadjust(mxAll / max(mxAll, [], 'all'), mice_contrast_2D{pts})); title(strcat(Exp_labels{m_num}{pts}, ' bsline-Max'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'fontsize', 14); 
    
    subplot(133); imagesc(imadjust(stdAll / max(stdAll, [], 'all'), mice_contrast_2D{pts})); title(strcat(Exp_labels{m_num}{pts}, ' bsline-Std'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'fontsize', 14); 
    
    set(gcf, 'Position',  [50, 50, 1500, 320]);
        nextPat = strcat('MxMnStd_VascularMap', Exp_labels{m_num}{pts}, '_p', num2str(m_num), '.png');
        saveas(gcf, nextPat);
    
end


%% Plot global raw and filtered time series
close all;

M_num = 1;
nextMse_data_raw = miceData_all{M_num}{1};
nextMse_data_mc = miceData_all{M_num}{2};
next_times = 1 : size(nextMse_data_mc{1}, 3);

bsln_color = [0 0 0.1];
stim_color = [1 0 1];

for pts = 1 : size(nextMse_data_mc,2)
    
    raw3D = nextMse_data_raw{pts};
    mc3D = nextMse_data_mc{pts};

    bslneRng = baseline_Ranges_all{pts};
    glo_raw_pct = getPercentChange_1D(squeeze(mean(mean(raw3D))), bslneRng);
    glo_mc_pct = getPercentChange_1D(squeeze(mean(mean(mc3D))), bslneRng);

    figure;
   
    set(gca,'fontsize',15); 
    title(strcat('Global Mean - ', Exp_labels{M_num}{pts}));
%     
%         patch([baseline_Ranges_all{pts}(1) baseline_Ranges_all{pts}(end) baseline_Ranges_all{pts}(end) baseline_Ranges_all{pts}(1)]/60, ...
%         [-100 -100 100 100], bsln_color, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
%     
%         patch([Stim_ranges{1}(end) Stim_ranges{6}(end) Stim_ranges{6}(end) Stim_ranges{1}(end)]/60, ...
%         [-100 -100 100 100], stim_color, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    
    hold on;
    plot(next_times/60, glo_raw_pct, '-b');
    plot(next_times/60, glo_mc_pct, '-r');

    grid;
    xline((baseline_Ranges_all{pts}(1))/60,'--k',{'baseline'}, 'fontsize',15); 
    xline(baseline_Ranges_all{pts}(end)/60,'--r',{'MK-injection'}, 'fontsize',15);
    ylabel('% pD signal change');
    xlabel('t (min)');
    yline(0,'-k',{''}); 
    ylim([-20 20]);
    xlim([3, 60]);
    legend( 'meanGlobal-Raw','meanGlobal-Filtered', 'Location','best');
    hold off;
    
    set(gcf, 'Position',  [50, 50, 900, 550]);
    
end

%% Auxiliary functions 

function dataCells_filtered = preProcess_mkDoseData(dataCell, params, strt, isTempRng)
    
    dataCells_filtered = {};
    for i = strt : size(dataCell, 2)
        
        miceData = {};    
        mc_3D_files = {}; raw_3D_files = {};
        nextPatient = dataCell{i};
        
        for a = 1 : size(nextPatient, 2)
            
            nextFN = nextPatient{a};
            nextData = h5read(nextFN, '/Data'); 
            next_times = h5read(nextFN, '/acqMetaData/time'); 
            
            nextData_raw_perm = permute(squeeze(nextData), [2 1 3]);
            
            if ~exist('isTemp','var')
                nextData_mc = fUSi_animalData_preProcessing(nextData, params); 

                nextScanFN = strcat(nextFN(1 : end - 5), '_HF_mc.scan');
                h5write(nextScanFN, '/Data', double(nextData_mc)); 
            else
                bslnRng = isTempRng;
                baselineTemplate = mean(nextData_raw_perm(:, :, bslnRng), 3);
                nextData_mc = fUSi_animalData_preProcessing(nextData, params, baselineTemplate); 

                nextScanFN = strcat(nextFN(1 : end - 5), '_HF_mc.scan');
                h5write(nextScanFN, '/Data', double(nextData_mc)); 
            end
                nextData_mc_perm = permute(squeeze(nextData_mc), [2 1 3]);

            raw_3D_files{a} = double(nextData_raw_perm);
            mc_3D_files{a} = double(nextData_mc_perm);

        end
        
        miceData{1} = raw_3D_files; 
        miceData{2} = mc_3D_files; 

        dataCells_filtered{i} = miceData;
        
    end
 end
 

 function pct = getPercentChange_1D(dat, brng)

    bsln = mean(dat(brng));
    ypct = 100*((dat - bsln) / bsln);
    
    pct = ypct;
 end

