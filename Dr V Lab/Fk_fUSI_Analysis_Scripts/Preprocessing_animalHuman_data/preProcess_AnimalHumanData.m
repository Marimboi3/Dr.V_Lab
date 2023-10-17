% ***** Function takes filenames of raw .scan data in a cell  and returns raw, filtered, 3D data and times
% after rigid, drift, or/and low pass filtering
% dataCell = raw .scan filenames
% params = [rgd p_rgd dtrn p_dtrn brd p_brd medn p_medn]


function dataCells_filtered = preProcess_AnimalHumanData(dataCell, params, strt, isTempRng)
    
    dataCells_filtered = {};
    for i = strt : size(dataCell, 2)
        
        patientData = {};    
        mc_3D_files = {}; raw_3D_files = {}; scs_times = {};
        nextPatient = dataCell{i};
        
        for a = 1 : size(nextPatient, 2)
            
            nextFN = nextPatient{a};

            if(nextFN(end-3:end)== 'scan')

                nextData = h5read(nextFN, '/Data'); 
                next_times = h5read(nextFN, '/acqMetaData/time'); 
            end

            if(nextFN(end-3:end)== '.acq')

                Fn_load = load(nextFN, '-mat');
                nextData = Fn_load.Acquisition.Data; 
                next_times = Fn_load.Acquisition.T;
            end
            
            nextData_raw_perm = permute(squeeze(nextData), [2 1 3]);
            
            if ~exist('isTemp','var')
                nextData_mc = fUSi_HumanData_preProcessing(nextData, params); 
            else
                bslnRng = isTempRng;
                baselineTemplate = mean(nextData_raw_perm(:, :, bslnRng), 3);
                nextData_mc = fUSi_HumanData_preProcessing(nextData, params, baselineTemplate); 
            end
                nextData_mc_perm = permute(squeeze(nextData_mc), [2 1 3]);


            raw_3D_files{a} = double(nextData_raw_perm);
            mc_3D_files{a} = double(nextData_mc_perm);
            scs_times{a} = double(next_times);

        end
        
        patientData{1} = raw_3D_files; 
        patientData{2} = mc_3D_files; 
        patientData{3} = scs_times; 

        dataCells_filtered{i} = patientData;
        
    end            
 end
