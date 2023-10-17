%**** Function save 2D image sequence into a movie  - animal specific ****
function fUSimageToMovie_enhanced(fname, dta, slen, map, mTitle, sFlag, contrast)

    if(size(dta, 2)  == 1)
        fData = permute(squeeze(dta), [2 1 3]);
    end
    if(size(dta, 2) > 1)
        fData = dta;
    end
    nmx = max(max(max(fData)));
    fData_sc = fData/nmx;
    
    figure;
    colormap(map);
   for si = 1 : slen
         slice_adjt = imadjust(fData_sc(:, :, si), contrast);
         imagesc(slice_adjt);
         title(mTitle);
         Mv(si) = getframe(gcf);
%          pause(0.1);
   end
   
   if (sFlag == 1)
       % StimON = VideoWriter(fname, 'Uncompressed AVI');
        StimON = VideoWriter(fname);
        StimON.FrameRate = 20;
        open(StimON);
        writeVideo(StimON, Mv);
        close(StimON);
   end
   
end