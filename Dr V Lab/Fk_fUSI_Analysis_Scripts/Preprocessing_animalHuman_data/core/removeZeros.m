function [dataOut, zeroRows, zeroCols] = removeZeros(dataIn, zeroRowsIn, zeroColsIn)

% This function removes any rows and any columns from fUS data that has
% zeros at all locations within the row/column. Alternativiely, if vargin
% are manually passed for zeroRows and zeroCols, those indices are removed
% instead.
%
% This function's original intended purpose was to remove any columns or
% rows that experienced clipping due to registration/motion correction of
% fUS data.
%
% function prototype:
% [dataOut, zeroRows, zeroCols] = removeZeros(dataIn, zeroRows, zeroCols)

[nRows, nCols, epochLength, nTrials] = size(dataIn); 

%% initialize
zeroRows = zeros(nRows,1);
zeroCols = zeros(1,nCols);

% for every frame, tally up any row/col that has a zero
for epoch = 1:epochLength
    for trial = 1:nTrials
        % this is the current frame
        tmpData = squeeze(dataIn(:,:,epoch,trial));
        % get indices for any row/col with all zeros
        tmpZeroRows = all(~tmpData,2);
        tmpZeroCols = all(~tmpData,1);
        % add any new indices to the overall indices
        zeroRows = any([zeroRows tmpZeroRows],2);
        zeroCols = any([zeroCols;tmpZeroCols],1);
    end
end

%% create data out by deleting zero rows/cols
dataOut = dataIn;

% find zero rows & columns (if not given)
if exist('zeroRowsIn','var') && exist('zeroColsIn','var')
    % Remove zero rows
    dataOut( zeroRowsIn , : , : , : ) = [];
    % Remove zero columns
    dataOut( : , zeroColsIn , : , : ) = [];
else    
    % Remove zero rows
    dataOut( zeroRows , : , : , : ) = [];
    % Remove zero columns
    dataOut( : , zeroCols , : , : ) = [];
end

return