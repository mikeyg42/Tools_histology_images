function [FINmasks, rectums, allROIsamples] = choosingROIs(outOfBoundsMask, imgCurrent, varargin)



%% ensure all your information is updated and available
close all force

if isstruct(varargin{1})
    handy = varargin{1};
else 
    handy = struct('ROISIZE', [], 'NUMROI', [], 'scaleFactor', [], 'directories', []);
    
    %how big would you like each ROI to be?
    handy.ROISIZE = [700, 700]; %this is in real world pixel units (same as RAW images, before resampling might've taken place
    
    %how many ROI's would you like to segment per region?
    handy.NUMROI = 10; 

    handy.scaleFactor = 1.75;  
    
    handy.rawdata_format = {'.tiff','.tif'};

    directoryImages = '/Users/jglendin/Desktop/tiled images as tiffs/';
    directorytosaveinto = '/Users/jglendin/Dropbox - Michael/Dropbox/processedMSimages/';
    
    directories = struct('saveDir', directorytosaveinto, 'imageDir', directoryImages);
    handy.directories = directories;
    
    imageDS = imagesdatastore(directoryImages, 'FileExtensions', rawdata_format);
end
    
 %% if handy is not supplied, it is assumed that this function is being called standalone, and not embedded in script
    %therefore we will not assume images are preprocessed!
    if ~isfloat(imgCurrent)
        imgCurrent = ensureDoubleScaled(imgCurrent);
% NOTE: if imgCurrent was saved from my imSegmentation_mosaics.m
% script, you should double check if it has already been downsampled!!!
        imgCurrent = preprocessRGBim(imgCurrent, scaleFactor);
    end
    
    if ~islogical(outOfBoundsMask)
        outOfBoundsMask > 0.5 = 1;
        outOfBoundsMask < 0.5 = 0;
        outOfBoundsMask = logical(outOfBoundsMask);
        outOfBoundsMask = imresize(outOfBoundsMask, size(imgCurrent, 1:2));
    end


%%



%<not done>




figyah = figure; 
set(figyah,'doublebuffer','off');

imDisp = imshow(imgCurrent);
set(gca, 'xlimmode','manual',...
'ylimmode','manual',...
'zlimmode','manual',...
'climmode','manual',...
'alimmode','manual',...
'Units', 'pixels');

%% start of loop to generate 10 ROI's

allROIsamples= zeros(size(imgCurrent, 1:2));
numSubsets=1;

while numSubsets ~= (1+NUMROI)

    win = randomCropWindow2d(size(imgCurrent), ROISIZE);
    rows = win.YLimits(1):win.YLimits(2);
    columns = win.XLimits(1):win.XLimits(2);
    rectangleROI = images.roi.Rectangle(imDisp.Parent, 'Position',[columns(1) rows(1) ROISIZE(2) ROISIZE(1)], 'Color', 'r');
    maskRect = createMask(rectangleROI);
    
    if sum(maskRect(:)) == 0
        delete(rectangleROI);
        continue
    end
    
    MASKintersect = maskRect.*outOfBoundsMask;
    overlappingROIs = maskRect.*allROIsamples;
    % here i say only keep the random window (ie "win") if it does not intersect with
    % the outofBoundsmask Or the already been sampled mask
    
    if sum(MASKintersect(:)) ~= sum(maskRect(:))
        delete(rectangleROI);
        continue
    end
    
    if sum(overlappingROIs(:)) > 0
        delete(rectangleROI);
        continue
    end
    
    allROIsamples = double(allROIsamples)+double(maskRect);
    FINmasks{numSubsets} = maskRect;
    rectums(numSubsets).points = rectangleROI.Vertices;
    
    numSubsets = numSubsets+1;
    
    delete(rectangleROI);
    
end