function segment_SaveFunc(mySettings, info, finalbinaryMask, imAdjRGB)
% segment_SaveFunction(mySettings, info, finalbinaryMask, imAdjRGB)
% the input "info" should be the imfinfo output when called on the original raw data; we
% only need this to get the dimensions of the raw data. Taking the reciprocal of the scale
% factor will introduce rounding errors that make your resized final images off by a pixel.

% Michael Glendinning, 2023

%/^*^\<--%--->\v.v/<--%->/^*^\<--%--->\v.v/<--%->/^*^\<--%--->\v.v/<--%->/^*^\<--%
%--%--\v.v/<--%->/^*^\<--%--->\v.v/<--%->/^*^\<--%--->\v.v/<--%->/^*^\<--%--->\v.v/

    Filename = split(mySettings.activeFilename, '/');
    activeFilename_wExt = split(Filename{end}, '.');
    activeFilename = activeFilename_wExt{1};


    mask_FMT = mySettings.fileFormats.segm_saveFMT_adjImage;
    RGB_FMT = mySettings.fileFormats.segm_saveFMT_foregroundMask;
    
    scaleFactor = mySettings.seg.seg_scaleFactor;
    USETHIS_scaleFactor = 1/scaleFactor;

    if mySettings.seg.seg_resizeBig
        finalSize = [info.Height, info.Width];
        finalbinaryMask = imresize(finalbinaryMask, finalSize, {@oscResampling, 4});
        imAdjRGB = imresize(imAdjRGB, finalSize, {@oscResampling, 4});
        
        try
            imwrite(finalbinaryMask, fullfile(mySettings.directories.saveSegm_foregroundMask,strcat(activeFilename,'scale_1.0_', '_ForegroundSegmented',mySettings.savePrefs.seg_versionID ,mask_FMT)), extractAfter(mask_FMT, '.'),...
                'Compression','none');
            imwrite(imAdjRGB, fullfile(mySettings.directories.saveSegm_adjImages, strcat(activeFilename,'scale_1.0_', '_adjustedRGBImage', mySettings.savePrefs.seg_versionID ,RGB_FMT)), extractAfter(RGB_FMT, '.'), ...
                'Compression','none');
        catch
            imwrite(imresize(finalbinaryMask,0.8, {@oscResampling, 4}) , fullfile(mySettings.directories.saveSegm_foregroundMask,strcat(activeFilename,'scale_0.8_', '_ForegroundSegmented',mySettings.savePrefs.seg_versionID , mask_FMT)), extractAfter(mask_FMT, '.'),...
                'Compression','none');
            imwrite(imresize(imAdjRGB, 0.8, {@oscResampling, 4}), fullfile(mySettings.directories.saveSegm_adjImages, strcat(activeFilename,'scale_0.8_', '_adjustedRGBImage', mySettings.savePrefs.seg_versionID, RGB_FMT)), extractAfter(RGB_FMT, '.'), ...
                'Compression','none');
        end
        
    else
        imwrite(finalbinaryMask, fullfile(mySettings.directories.saveSegm_foregroundMask,strcat( activeFilename,'scale_', sprintf('%.4g', USETHIS_scaleFactor), '_ForegroundSegmented',mySettings.savePrefs.seg_versionID,mask_FMT)), extractAfter(mask_FMT, '.'),...
            'Compression','none');
        
        imwrite(imAdjRGB, fullfile(mySettings.directories.saveSegm_adjImages, strcat( activeFilename,'scale_', sprintf('%.4g', USETHIS_scaleFactor), '_adjustedRGBImage',mySettings.savePrefs.seg_versionID, RGB_FMT)), extractAfter(RGB_FMT, '.'), ...
            'Compression','none');
    end
end