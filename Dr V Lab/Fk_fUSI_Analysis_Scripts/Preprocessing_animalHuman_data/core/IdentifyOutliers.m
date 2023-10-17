function OutlierIndex = IdentifyOutliers(data,method,varargin)

if any(strcmp(varargin,'verbose'))
    verbosity = varargin{find(strcmp(varargin,'verbose'))+1};
else
    verbosity = false;
end

if any(strcmp(varargin,'threshold'))
    Threshold = varargin{find(strcmp(varargin,'threshold'))+1};
else
    switch method
        case 'automatic'
            Threshold = 5;
        case 'manual'
            Threshold = 100;
    end
end

data_resized=reshape(data,1,[]);


switch method
    case 'automatic'
        OutlierIndex = isoutlier(data_resized,'median','ThresholdFactor',Threshold);
    case 'manual'
        OutlierIndex = abs(data_resized)>=Threshold;
end

OutlierIndex = reshape(OutlierIndex,size(data));

if verbosity
    if any(OutlierIndex)
        fmt = ['outliers removed: ', repmat('%g, ', 1, numel(data(OutlierIndex))-1), '%g \n'];
        fprintf(fmt, data(OutlierIndex));
    end
end

end %End of function