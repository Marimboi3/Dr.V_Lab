% uninterp saved data is a script that will run through a chosen directory,
% prompt the user to choose which files need uninterpolating, load those
% files, uninterp the doppler data, and overwrite the existing files. This
% script was created in November 2020 as part of the new core that does not
% use the interp2(data,2) call in create doppler images. This is a faster
% method of converting to the new core method than actually re-running the
% core, which requires re-beamforming all those data.

%% get the sessions & runs desired

% get the directory wanted
folder_path = uigetdir(getUserDataPath('pathType','processed'));
if ~folder_path, return, end

% get the files in that directory
files = struct2table(dir(folder_path));
doppler_filenames = files.name(contains(files.name,'doppler'));

% use a GUI to select which ones you want
[indx,tf] = listdlg('PromptString','Select Session/Run(s) to load:',...
    'SelectionMode','multiple','ListString',doppler_filenames);
update_list = doppler_filenames(indx);

%% update files
for i = 1:length(update_list)
    % get the filename
    filename = fullfile(folder_path,update_list{i});
    fprintf('Updating file %i/%i (%s)\n',i,length(update_list),filename)  
    try
        % load the file
        fprintf('Loading...'), tmp = load(filename); fprintf('done.')
        % check to see if it's already the right size (skip if so)
        if size(tmp.iDop,1) <= 128
            fprintf('This Doppler data is already uninterpolated. Skipping. (size: %i x %i)\n',...
                size(tmp.iDop,1), size(tmp.iDop,2))
        else
            % uninterpolate
            fprintf('\tUninterpolating...')
            tmp.iDop = preProcess(tmp.iDop,'uninterp',true);
            tmp.wallpaper = preProcess(tmp.wallpaper,'uninterp',true);        
            % save
            fprintf('done (new size: %i x %i).\tSaving...',size(tmp.iDop,1),size(tmp.iDop,2))
            save(filename,'-struct','tmp')   
            fprintf('done.\n')
        end        
    catch
        warning('%s failed.',filename)
    end
end