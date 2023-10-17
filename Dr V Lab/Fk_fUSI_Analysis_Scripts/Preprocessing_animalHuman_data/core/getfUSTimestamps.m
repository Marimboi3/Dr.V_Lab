function timestamps = getfUSTimestamps(Session, Run)


%% find fUS acqs times (Verasonics' log file)
cd(getUserDataPath('pathType','fUS','session',Session,'run',Run));
t_fUSacqs = abs(textread('timeLog.txt', '%f')); %#ok<DTXTRD> % load fUS acq time log

% reshape fUS time log to [HH MM SS.MS]
if Session < 30 || Session == 35 || (Session == 39 && Run ~= 1)
    % older style fUS time log has 4 columns
    t_fUSacqs = reshape(t_fUSacqs,4, [])';  
    t_fUSacqs = t_fUSacqs(:,1:3);                   
else
    % new fUS time log has 7 columns, 4-6 are already sync'd to behavior ctrlr
    t_fUSacqs = reshape(t_fUSacqs,7, [])';
    t_fUSacqs = t_fUSacqs(:,4:6);                   
end

timestamps = t_fUSacqs(:,1)*3600 + t_fUSacqs(:,2)*60+ t_fUSacqs(:,3); % convert to [s]
end