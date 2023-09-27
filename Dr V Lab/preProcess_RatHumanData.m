% ***** Function takes filenames of raw .scan data in a cell  and returns raw, filtered, 3D data and times
% after rigid, drift, or/and low pass filtering
% dataCell = raw .scan filenames
% params = [rgd p_rgd dtrn p_dtrn brd p_brd medn p_medn]

function [raw3D, mc3D, scsTimes] = preProcess_RatHumanData(dataCell, params)
    mc_3D_files = {}; raw_3D_files = {}; scs_times = {};
    
    for a = 1 : size(dataCell, 2)
        
        nextFN = dataCell{a};
        nextData = h5read(nextFN, '/Data'); 
        next_times = h5read(nextFN, '/acqMetaData/time'); 
        
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