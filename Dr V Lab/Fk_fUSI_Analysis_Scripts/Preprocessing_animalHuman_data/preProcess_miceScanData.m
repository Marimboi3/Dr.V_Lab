
function dataCells_filtered = preProcess_miceScanData(dataCell, params, strt, isTempRng)
    
    dataCells_filtered = {};
    for i = strt : size(dataCell, 2)
        
        miceData = {};    
        mc_3D_files = {}; raw_3D_files = {};
        nextPatient = dataCell{i};
        
        for a = 1 : size(nextPatient, 2)
            
            nextFN = nextPatient{a};
            nextData = h5read(nextFN, '/Data'); 
            next_times = h5read(nextFN, '/acqMetaData/time'); 
            
            nextData_raw_perm = permute(squeeze(nextData), [2 1 3]);
            
            if ~exist('isTemp','var')
                nextData_mc = fUSi_animalData_preProcessing(nextData, params); 

                nextScanFN = strcat(nextFN(1 : end - 5), '_mc.scan');
                h5write(nextScanFN, '/Data', double(nextData_mc)); 
            else
                bslnRng = isTempRng;
                baselineTemplate = mean(nextData_raw_perm(:, :, bslnRng), 3);
                nextData_mc = fUSi_animalData_pProcessing(nextData, params, baselineTemplate); 

                nextScanFN = strcat(nextFN(1 : end - 5), '_mc.scan');
                h5write(nextScanFN, '/Data', double(nextData_mc)); 
            end
                nextData_mc_perm = permute(squeeze(nextData_mc), [2 1 3]);

            raw_3D_files{a} = double(nextData_raw_perm);
            mc_3D_files{a} = double(nextData_mc_perm);

        end
        
        miceData{1} = raw_3D_files; 
        miceData{2} = mc_3D_files; 

        dataCells_filtered{i} = miceData;
        
    end
 end
 