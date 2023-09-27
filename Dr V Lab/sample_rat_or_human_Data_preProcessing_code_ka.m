% --------- Sample fUSI data pre-processing pipeline - Rat or Human data -----------------

%% Initialize workspace
clear all; close all; 
clearvars;
clc;

%% Init fUS SC data .scan file names and labeals

Fname_1 = 'sample_data_1.scan';      % File names - can be any number of files (Files must be located in a folder in the the matlab path)
Fname_2 = 'sample_data_2.scan';
Fname_3 = 'sample_data_3.scan';
Fname_4 = 'sample_data_4.scan';
Fname_5 = 'sample_data_5.scan';

FN_all_cell = {Fname_1 Fname_2 Fname_3 Fname_4 Fname_5};         % Create cell with file names
FN_labels = {' FLabel-1' ' FLabel-2' ' FLabel-3' ' FLabel-4' ' FLabel-5'};                                                           % Create cell of file labels

%% Filter/Pre-process all patient fUSI data 
sTart = 1;
lp1 = 0.0001;
lp2 = 0.05;     %upper limit of the low pass filter - can be adjusted
Fq = 1;

parameters = [1, 0, 0, 40, 1, [lp1, lp2, Fq], 0, 11];   % params = [rgd p_rgd dtrn p_dtrn lwps p_lwps medn p_medn]
[dataCell_3D_raw, dataCell_3D_mc, dataCell_times] = preProcess_RatHumanData(FN_all_cell, parameters, sTart);

% Save data cells if needed - You might need to create a dummy file with the file name before hand
% save('sampleExperiment_data_rawFiltered.mat', 'dataCell_3D_raw', 'dataCell_3D_mc', 'dataCell_times','-append');

%% Visualize raw and filteres fUSI SC Data
contrast_raw = {[0.0 0.1] [0.0 0.1]  [0.0 0.1]  [0.0 0.1]  [0.0 0.1] };
contrast_mc = {[0.0 0.1] [0.0 0.1]  [0.0 0.1]  [0.0 0.1]  [0.0 0.1] };

%% Show short movie of the raw and filtered data
close all;
FN_num = 1;  % patient number
mov_len = 100; % movie length (sec) to display
next_patient = FN_labels{FN_num};

% *** View raw data
title_raw = strcat(next_patient, '  Raw'); 
raw_FN = strcat(next_patient, '_raw'); 
raw_data = dataCell_3D_raw{FN_num};
fUSimageToMovie_enhanced(raw_FN, raw_data, mov_len, 'hot', title_raw, 0, contrast_raw{FN_num});  % 1=save movie, 0 = don't save 

% *** View filtered data
title_filt = strcat(next_patient, '  Filtered'); 
filt_FN = strcat(next_patient, '_mc'); 
filtered_data = dataCell_3D_mc{FN_num};
fUSimageToMovie_enhanced(filt_FN, filtered_data, mov_len, 'hot', title_filt, 0, contrast_mc{FN_num});  % 1=save movie, 0 = don't save 

%% Show 2D power Doppler mean/max  figures (Can save file)
close all;
contrast_2D = {[0.0 0.2] [0.0 0.2]  [0.0 0.2]  [0.0 0.2]  [0.0 0.2] }; % contrast level for each data

for pts = 1 : size(dataCell_3D_mc,2)
    
    mc3D = dataCell_3D_mc{pts};                       % next mc data
    dataRng = 5 : size(mc3D, 3) - 5;                    % data range
    
    mxAll = max(mc3D(:, :, dataRng), [], 3);        % 2D max in time dim
    mnAll = mean(mc3D, 3);                               % 2D minimum
    
    figure;
    colormap(hot);
    subplot(121); imagesc(imadjust(mnAll / max(mnAll, [], 'all'), contrast_2D{pts})); title(strcat(FN_labels{pts}, ' Mean all'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'fontsize', 14); 
    subplot(122); imagesc(imadjust(mxAll / max(mxAll, [], 'all'), contrast_2D{pts})); title(strcat(FN_labels{pts}, ' Max all'));
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'fontsize', 14); 
    set(gcf, 'Position',  [50, 50, 1220, 400]);
    
    % nextPat = strcat('epid_MxMn', FN_labels{P_num}, '.png');  % Can save images of desired
    % saveas(gcf, nextPat); 
    
end

%% Data analysis
acqPeriod = 1;                                                   % fUSI acquisition period (s)
acqFreq = 1 / acqPeriod;   
baseline_Ranges = {1:30 1:30 1:30 1:30 1:30};            % specify baseline range for each data

%% Plot global raw and filtered time series
close all;
stim_color = [1 0 1];

for pts = 1 : size(dataCell_3D_raw,2)
    
    raw3D = dataCell_3D_raw{pts};
    mc3D = dataCell_3D_mc{pts};
    next_times = dataCell_times{pts};

    bslneRng = baseline_Ranges{pts};
    glo_raw_pct = getPercentChange_1D(squeeze(mean(mean(raw3D))), bslneRng);
    glo_mc_pct = getPercentChange_1D(squeeze(mean(mean(mc3D))), bslneRng);

    figure;
   
    set(gca,'fontsize',15); 
    title(strcat('Global Mean - ',FN_labels{pts}));
    plot(next_times/60, glo_raw_pct, '-b');
    hold on;
    plot(next_times/60, glo_mc_pct, '-r');

    grid;
    ylabel('% pD signal change');
    xlabel('t (min)');
    yline(0,'-k',{''}); 
    ylim([-100 100]);
    legend('meanGlobal-Raw','meanGlobal-Filtered', 'Location','south');
    hold off;
    
end

%% auxiliary function

 function pct = getPercentChange_1D(dat, brng)

    bsln = mean(dat(brng));
    ypct = 100*((dat - bsln) / bsln);
    
    pct = ypct;
end


