%***** General Function for Image and Results Plots *************

% Main
function VisualizePlot(ptype, xdta, ydta, params)
    if(ptype == 'slices')
        fus_image_slices(ydta);
    end
end

% *** Visualize 2D fUS time series 
function fus_image_slices(img)
     f_bline = figure;
    for ai = 1 : size(img, 3)
        Img_slice = img(:, :, ai);
        figure(f_bline); 
        imagesc(Img_slice);
        colormap(turbo);
%         pause(0.1);
        drawnow;
    end
end