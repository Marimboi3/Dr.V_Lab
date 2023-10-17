%
% **** Function save 2D image sequence into a movie  ****

function imageSeqToMovie(fname, dta, slen, map, mTitle, sFlag)
    fData = permute(squeeze(dta), [2 1 3]);
    
    figure;
    colormap(map);
   for si = 1 : slen
         imagesc(fData(:, :, si));
         title(mTitle);
         Mv(si) = getframe(gcf);
   end
   
   if (sFlag == 1)
       % StimON = VideoWriter(fname, 'Uncompressed AVI');
        StimON = VideoWriter(fname);
        StimON.FrameRate = 30;
        open(StimON);
        writeVideo(StimON, Mv);
        close(StimON);
   end
   
end
