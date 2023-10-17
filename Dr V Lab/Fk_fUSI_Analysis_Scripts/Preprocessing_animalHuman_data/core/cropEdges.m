function dopOut = cropEdges(dopIn, edge)
% function dopIn = cropEdges(dopIn, edge)
%
% varin:
% dopIn is a [yPixels, xPixels, nWindows, nTrials] array of doppler data
% typically loaded from the core processing data files, e.g. using
% loadSingleRunData or manually loading in doppler_S##R##.mat.
%
% edge is the number of pixels you would like to shave from each edge
% For example, setting edge to 50 would shave 50 pixels from every edge in
% the image. In a Y x X image, result is (edge+1:Y-edge) x (edge+1:X-edge).

%% print data in
[yPix, xPix, ~, ~] = size(dopIn);
fprintf('Doppler data contains %i x %i pixels. Removing %i pixels from edges... ', yPix, xPix, edge);

%% check for errors
if ((yPix - 2*edge)<2) || ((xPix - 2*edge)<2)
    error('The image is not large enough to remove this many voxels.')
end

%% perform operation
yPixelsOut = edge+1 : yPix-edge;
xPixelsOut = edge+1 : xPix-edge;
dopOut = dopIn(yPixelsOut, xPixelsOut,:,:);

%% print data out
[yPixOut, xPixOut, ~, ~] = size(dopOut);
fprintf('Doppler out is %i x %i.\n', yPixOut, xPixOut);