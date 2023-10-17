function numTime = parseTime(chartime,columnTitle) %#ok<INUSL>
%% Converts the character arrays in the behavior data to numbers
% Moved from being a subfunction in parseBehavior to being an independent
% function due a restructuring of parseBehavior requiring this function in
% multiple different functions
%
% Written by Sumner (I think)

chartime = eval(['cell2mat({chartime.' columnTitle '}'');']);

% extract hours/min/sec/ms from chartimes
%#ok<*ST2NM>
HH = str2num(chartime(:,1:2));   % Hours array
MM = str2num(chartime(:,4:5));   % Minutes array
SS = str2num(chartime(:,7:8));   % seconds array
MS = str2num(chartime(:,10:12)); % milliseconds array

% concat visual stim trial times [HH MM SS.MS]
numTimeArray = cat(2, HH, MM, SS+1e-3*MS);

% convert time to [s]
numTime =   numTimeArray(:,1)*3600 + numTimeArray(:,2)*60+ numTimeArray(:,3);
end