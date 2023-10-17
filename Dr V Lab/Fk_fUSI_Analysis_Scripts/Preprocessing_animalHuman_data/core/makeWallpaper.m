function wallpaper = makeWallpaper(DopplerImage)
    wallpaper = nthroot(squeeze(mean(mean(DopplerImage,3),4)),3);
    lowerCutoff = quantile(wallpaper(:),0.01);
    wallpaper(wallpaper<lowerCutoff) = lowerCutoff;
end