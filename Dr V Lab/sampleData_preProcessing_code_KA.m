%
% --------- Sample pre-processing and data analysis code - Urodynamics  (Standard Protocols) -----------------
%

%% Initialize and Load Files
clear all; close all; 
clearvars;
clc;

%% load Rat .scan data

Fn_Rat7_1 = 'R7_0_1_Sag.scan';
ratData = h5read(Fn_Rat7_1, '/Data'); 

%% Pre-process to remove motion artifacts
lp1 = 0.0001;   % start frequency
lp2 = 0.02;      % End frequency
Fq = 1;
parameters = [1, 0, 0, 40, 1, [lp1, lp2, Fq], 0, 11];   % params = [rgd p_rgd dtrn p_dtrn lwps p_lwps medn p_medn]

ratData_mc = fUSdata_psrd_filter(ratData, parameters , 0); % Last input: 1 = perform drift correction first, 0 = after rgid motion correction 
ratData_raw_3D = permute(squeeze(ratData), [2 1 3]);
ratData_mc_3D = permute(squeeze(ratData_mc), [2 1 3]);

%% View movie of raw and filtered 3D data

close all;
mov_len = 150; 
contrast = [0 0.03];
next_animal = 'sample-rat';
% *** View raw data
title_raw = strcat(next_animal, '  Raw'); 
raw_FN = strcat(next_animal, '_raw'); 
raw_data = ratData_raw_3D;
fUSimageToMovie_enhanced(raw_FN, raw_data, mov_len, 'hot', title_raw, 0, contrast);  % 1=save movie, 0 = don't save 

% *** View filtered data
title_filt = strcat(next_animal, '  Filtered'); 
filt_FN = strcat(next_animal, '_mc'); 
filtered_data = ratData_mc_3D;
fUSimageToMovie_enhanced(filt_FN, filtered_data, mov_len, 'hot', title_filt, 0, contrast);  % 1=save movie, 0 = don't save 

%% 2D mean/max power doppler figures
close all
Rat_contrast = [0.0 0.07];

mc3D = ratData_mc_3D;
dataRng = 100 : size(mc3D, 3) - 100;
slc = 371;

mxAll = max(mc3D(:, :, dataRng), [], 3);
mnAll = mean(mc3D, 3);

figure;
set(gca,'fontsize', 15); 
colormap(hot);
subplot(121); imagesc(imadjust(mnAll / max(mnAll, [], 'all'), Rat_contrast)); title(strcat('SampleRat', ' Mean all'));
set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
subplot(122); imagesc(imadjust(mxAll / max(mxAll, [], 'all'), Rat_contrast)); title(strcat('SampleRat', ' Max all'));
set(gcf, 'Position',  [50, 50, 1200, 400]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);

%% Data analysis
