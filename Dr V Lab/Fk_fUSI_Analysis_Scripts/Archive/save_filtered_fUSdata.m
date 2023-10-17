% **** Function saves motion artifacts corrected fUS image data ****
% fName = .acq filename
% img4D = corrected fUS image data
% rwData = original data struct
function save_filtered_fUSdata(fName, img4D, rwData)

disp('Saving motion and breathing corrected fUS data....');
fUS_Perm = img4D;
% fUS_Perm = permute(img4D, [1 3 2 4]);
rwData.Acquisition.Data = fUS_Perm;
Acquisition = rwData.Acquisition;
save(fName, 'Acquisition');

disp('Done saving .acq file!');

end