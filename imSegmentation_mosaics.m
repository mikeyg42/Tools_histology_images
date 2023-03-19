function varargout = imSegmentation_mosaics(varargin)
% Syntax OPT 1: loopingSegmentation_mosaics (... no inut/output, instead modify the first few lines of this .m file)
% Syntax OPT 2: {mask1, mask2,..., maskn} = loopingSegmentation_mosaics {path1, path2,...pathn};
% This function is the parent function for all of the foreground/background 
% image segmentation code I've written for whole-slide image (WSI) mosaics.

% There are a few ways that it runs: 

% Opt 1. Standard use case is NO input and NO output variables.
% Instead, in the first chunk of the code in this .m file contains some file 
% paths that should be replaced to reflect your local storage. These will
% be your inputs and outputs effectively.
%      - Importantly, image filenames are assumed to contain substantial information
%          as to that image's etiology. It must be easily parsed and easily
%          locatable in the filesystem
%      - note: Its is a Terrrible idea to ever run code on your only copy of a dataset... 
%          Do yourself a favor, back it up, and then work on a copy.

% Opt 2. Atypical use case. INPUT is a cell array of file paths. OUTPUT will be a
% cell array of the same length as INPUT, except each cell will have a binary image files.  

% BACKGROUND -  This code was made to address new challenges associated with 
% my ever-increasing propensity of whole slide images for any chromogenic
% histology staining experiment. After ISH or IHC, I would image my tissue using 
% our brightfield scope (an Zeiss AxioImager), equipped w/ a 10X objective
% and automated stage, and Zeiss's proprietary software Zen Black 2012. I
% also used Zen for stitching of my images. 

% The functions I've included here SHOULD work just as well for HE, LFB,
% Nissl, ect. THE CRITICAL ASSSUMPTION is that the ENTIRE stained tissue
% is contained within the image, ie there should be a large blob in the
% middle of a white sea. I've built in some flexibility so that if your blob
% touches the edge once or twice its OKAY. But the majority of edges need to be
% clear for any of this to work. Some aspects of the code will also assume
% each WSI has 1 blob, but it is an ongoing effort to handle 2+ blobs. The 
% SECOND caveat is that will NOT segment fluroescent histolgoy.

%  The script, as is, relies heavily upon your file naming convention.
%  Prepending my convention with your own easy and will ensure you have no
%  issues. It requires three, underscore seperated identifiers: 

%               {1} _ {2} _ {3} _ {optional 4}.tiff
%        ie. "GROUPid_SAMPLEid_STAINid_ExtraInfo.tiff   

% {1} - experimental group ID. I often use: "SICK" or "HEALTHY". 
%       
% {2} - SAMPLE id. I often use A, B, C, ...

% {3} - identification of the particular staining performed in that image
%       (H+E, ect.) For IHC usually the particular target antigen for
%       IHC-Dab (e.g. CD31, Sox17, ...)

% {4} - whatever other notes/identifiers you want, it should not mess up
%       anything (...except potentially you shouldn't include more underscores
%       or repeat one of the identifiers from {1}, {2}, or {3}?)

% After the input is parsed, images are preprocessed. Raw image mosaics
% always exhibit a hallmark, frustrating grid-pattern (I'll refer to as
% the "picnic tablecloth effect"). The best way to prevent this is careful
% setting up of Kohler illumination and optimizing of image acquisition settings.
% But, after acquisition, your solution depends on how your downstream 
% analyses's sensitivity. Here I've addressed the picnic table effect with a 
% background correction algorithm (using the tophat filter). If this is 
% insufficient, you can modify preprocessingRawRGBims.m and backgroundCorrectRGB.m 
% And for more serious problems, check out the open source dockerized 
% workflow called MCMICRO, which uses the program BaSiC, (also an ImageJ plugin). 

% This script allows you to pick up where you left off.....




%<FINISH THIS>%



% Michael Glendinning, 2023

close force all
clear all

%% set your parameters

    %1. Set your scale factor. The closer to 1 the better! If you
    %experience crashing or painfully slow processing, raise this.
    scaleFactor = 2.8;  
    
    %2. Define the file formats in which raw images will use the following
    %extensions:1.4 or
    rawdata_format = {'.tiff','.tif'};

    %3. delete raw data files after processing and saving? 
    deleteYES = 1;
% NOTE: AFTER PROCESSING, ALL RAW DATA IN directoryImages WILL BE DESTROYED, iff deleteYES equals 1 **
    
    %4. set these directory locations for your local data storage
% NOTE: you MUST include the file seperator "/" at the end of each path...

    directoryImages = '/Users/jglendin/Desktop/tiled images as tiffs/';
    directorytosaveinto = '/Users/jglendin/Dropbox - Michael/Dropbox/processedMSimages/';
   
    %5. (optiona;) Include additional text to append to each save file (
    versionNum = '_v1';
    

    
%% OPTIONAL. start from somewhere in yoru directory of images other than the first file

counter = 2;

% note, counter allows us to iteratively navigate our image datastore with our data.
%       Images in the datastore are in the same order as they will be in
%       the field "filenames" within our master array "handy".(i.e.
%       handy.filenames{handy.counter} should always tell you which file is active.

%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%
%%  Don't change anything else 
%   |     |     |     |     |
%   V     V     V     V     V

%%  1. Load images and log some information about them
if nargin ~= 0 && nargout ~= 0
    if nargin ~= nargout
        disp('# of inputs must equal # of outputs! try again..');
        return
    end
    % create a datastore for all the images listed in input
    imageDS = imageDatastore(varargin);
    outputFlag = 1;
else 
    % OR create a datastore for all images found within image directory above
    imageDS = imageDatastore(directoryImages, 'FileExtensions', rawdata_format);
    outputFlag = 0;
end

% get the sample from the full filenames and load into a cell array - (called "imageNames")
fullFileNames = vertcat(imageDS.Files);
[~, imageNames, ~] = cellfun(@fileparts, fullFileNames, 'UniformOutput', false);

%% 2. Create a handy structure called "handy" 
% contains exp parameters, directory paths, tracks your progress, ect.

my_dirs = struct('loadDir', directoryImages, 'saveDir', directorytosaveinto);

handy = struct('counter', counter, 'versionNum', versionNum,...
    'scaleFactor', scaleFactor, 'images', imageDS, 'OutputFlag', outputFlag,...
    'deleteYES', deleteYES, 'directories', my_dirs);

handy.filenames = imageNames;

USETHISscaleFactor = 1/scaleFactor; 

savedfileList  = struct2cell(dir(handy.directories.saveDir))';
savedfileList = savedfileList(:,1); % isolate the "filename" column
allfiles_w_tiff = savedfileList(contains(savedfileList, 'tiff'));
alllfiles_w_png = savedfileList(contains(savedfileList, 'png')); % will be foreground/background segementation masks only
allSavedFiles = savedfileList(contains(savedfileList, versionNum )); % this gets rid of the hidden files named ".." and such

%--%--%--%--%--% Done w/ prep! --%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%

%%  3. load up 1st image and start processing

nImages = numel(imageNames);

while counter<nImages % loop to go through your entire image datastore
activeFilename = handy.filenames{handy.counter};

% prints the date and time into command window
printStartTime_imProcessing(activeFilename);

% I can leverage an image's format to inform me of what the contexts are :
%   RGB/grayscale images saved as ".tiffs",  binary segmentation saved as ".png"

% I have saved files in the format of: "ExpGROUPid_SAMPLEid_STAINid_anyExtraInfo.tiff" e.g. "MS_317_TMEM119_stitched19.tiff"
namePieces = split(activeFilename, '_');
groupID =  namePieces{1};
sampleID = namePieces{2};
stainID = namePieces{3};

% Now we look within the dir "allSavedFiles" for filenames sharing these ID's:
matchingIDandStain = contains(allSavedFiles, sampleID, 'IgnoreCase',true).*contains(allSavedFiles, stainID, 'IgnoreCase',true);

% How we then proceed depends on how many matches we find...
numberMatches = sum(matchingIDandStain);

switch num2str(numberMatches)
    
    % ...if we have zero matches, we start from beginning of processing
    % workflow
    case '0'
        % read in the raw image...
        [imgCurrentlyBig, ~] = readimage(handy.images, handy.counter); % handy.images is my raw image datastore!!!
        
        % preprocess RGB image: 1. apply the downsampling scaleFactor, 2. convert image to
        % float (double), 3. background/flatfield correction
        tic; imAdjRGB = preprocessRawRGBims(imgCurrentlyBig, USETHISscaleFactor);toc;
        
        % call the segmentation GUI function for foreground/background segmenting
        finalbinaryMask = segmentationWrapper(handy, imAdjRGB);
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        
        % ..if we have ONE match in save directory, we can skip ONE step!
        %   (but which step?...depends on if our match is the RGB or binary file)
    case '1'
        
        my1savedfile = allSavedFiles(find(matchingIDandStain));
        
        % OPTION 1- we only have saved a TIFF file
        if sum(contains(allfiles_w_tiff,my1savedfile))==1
            
            % use imfinfo function to ascertain raw image's dimensions 
            info = imfinfo(strcat(handy.directories.saveDir, my1savedfile));
            targetCols = info.Width*USETHISscaleFactor;
            targetRows = info.Height*USETHISscaleFactor;
            
            % tiff files = RGB, so read in the RGB image match from saveDir
            imAdjRGB = imread(strcat(handy.directories.saveDir, my1savedfile));
            imAdjRGB = imresize(imAdjRGB, [targetRows, targetCols, 3], {@oscResampling, 4});
            
            % call the segmentation GUI function for foreground/background segmenting
            finalbinaryMask = segmentationWrapper(handy, imAdjRGB);
            
            % OPTION 2 = we have saved a PNG file
        elseif sum(contains(alllfiles_w_png,my1savedfile))==1
            info2 = imfinfo(strcat(handy.directories.loadDir, directoryImages));

            sz = [info2.Height; info2.Width].*USETHISscalefactor;
            
            % png files = binary masks of foreground/background segmentation
            binaryMask = imread(strcat(handy.directories.saveDir, my1savedfile), 'png', 'BackgroundColor', 0);
            binaryMask = imresize(binaryMask, [sz(1), sz(2)], {@oscResampling, 4});
            binaryMask = logical(binaryMask);
            
            
            % read into workspace the raw image file, then the preprocessing
            % fcn to generate the adjusted RGBim. then we can proceed to
            % segmentation function
            [imgCurrentlyBig, ~] = readimage(handy.images, handy.counter);
            tic; imAdjRGB = preprocessRawRGBims(imgCurrentlyBig, USETHISscaleFactor); toc;
            finalbinaryMask = segmentationWrapper(handy, imAdjRGB, binaryMask);
            
            % OPTION 3 (unlikely) = we have some non-png, non-tiff fileformats in saveDIr
        else %ie you have one match in the directory, but its neither a tiff or a png.
            disp('you have errant fileformats in your saveDir. Get rid of anything not PNG and TIFF');
            return;
        end
        
        %-------------------------------------------------------------------------%
        %-------------------------------------------------------------------------%
        % 2 matches is the BEST case scenario! 2 matches suggests we have what we need:
        %     TWO images in our saveDir matching the identifiers in our active file.
        
    case '2'
        
        my2savedfiles = allSavedFiles(find(matchingIDandStain));
        
        % let's quickly double check it is 1x png and 1x tiff, and not 2x png or something....
        % SANITY CHECK:
        if sum(contains(alllfiles_w_png,my2savedfiles))~=1 || sum(contains(allfiles_w_tiff,my2savedfiles))~=1
            %if this is true, then save dir does not actually contain 1 image + 1 mask of your stain of interest
            disp(strcart('error!: your savedDir contains duplicates that are problematic... fix:', activeFilename))
            return
        end
        % end sanity check

        % read in the already adjusted RGB image
        tiffFile = strcat(handy.directories.saveDir, my2savedfiles{contains(my2savedfiles, 'tiff')==1});
        imAdjRGB = imread(tiffFile, 'tiff');
        
        % read in binary mask image
        pngFile = strcat(handy.directories.saveDir, my2savedfiles{contains(my2savedfiles, 'png')==1});
        binaryMask = imread(pngFile, 'png', 'BackgroundColor', 0);
        binaryMask = logical(binaryMask);
        
        % visualize + finalize foreground/background segmentation mask
        finalbinaryMask = segmentationWrapper(handy, imAdjRGB, binaryMask);
        
    otherwise
        disp(strcat('error!: >2 files in this dir match the active filename: ', activeFilename));
        return
        
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
        
%% 4. Saving!!




% CAUTION: this will overwrite whatever might have been there before!!
imwrite(finalbinaryMask, strcat(handy.directories.saveDir, activeFilename,'_ForegroundSegmented_v1.png'), 'png',...
    'Compression','none');

imwrite(imAdjRGB, strcat(handy.directories.saveDir, activeFilename,'_adjustedRGBImage_v1.tiff'), 'tiff', ...
    'Compression','none');


%% DONE! first save, then move on to the next image

% update the counter!
handy.counter = handy.counter+1;

end

end






function finalMask = segmentationWrapper(handy, imAdjRGB, varargin)
%% segmentation of the tissue piece from the background. 
% 
% Inputs    : 1. pre-processed, high resolution whole slide image
%             2. (optional) an already-made segmentation mask delineating foreground/background
%-------------
% Processes : 1. the formation of the initial tissue mask.
%             2. a la carte menu of methods to remake mask or refine one
%             3. saving the mask and the imAdjRGB locally for downstream
%-------------
% Outputs   : 1. segmented binary image (white=foreground)

% Highlights: 
% - a few  of the segmentation techniques are (as far as I can tell)
% entirely novel methodologies. 
% - The semi-automated approach is very powerful here, as one is forced to 
% become very well acquainted with both their data and their code, and 
% yet, if done well, minimal biases are introduced. 
%  - Some of these segmentation methods are not my own - I credit each
%  as applicable, although I've modified them all to some extent

% CAVEAT, for downstream applications involving high precision 
% quantification its inappropriate to pick and choose segmentation methods
% and refinements. Find what works best for your images and stick to that!

% Michael Glendinning, 2022

%% in case there is already a draft mask...
if numel(varargin) > 0 && ismatrix(varargin{1})
   binaryMask = varargin{1};
else
% ... otherwise segment the image initially using morphology
 binaryMask = morphologySegment(imAdjRGB);   
end
%% Mask-Refining GUI
finalMask = finalizeMASK(binaryMask,imAdjRGB);


end