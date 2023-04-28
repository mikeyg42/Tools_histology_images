function [mySettings] = setts_and_prefs

%% directories + file locations
%set these directory locations for your local data storage

% SEGMENTATION draws from "directOfImages" and saved saveDestination_adjImages 
% the corrresponding masks are saved in saveDestination_foregroundMask.

% REGISTRATION requires 2 images and 2 masks, which it gets from the 2 folders SEGMENT.
% saved into.

% ======================================================================================== %

% Directories!:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directOfImages          = '/Users/jglendin/Dropbox - Michael/Dropbox/';
saveSegm_adjImages      = '/Users/jglendin/Dropbox - Michael/Dropbox/processedMSimages/'; 
saveSegm_foregroundMask = '/Users/jglendin/Dropbox - Michael/Dropbox/processedMSimages/';
saveRegistraion         = '/Users/jglendin/Dropbox - Michael/Dropbox/processedMSimages/';
saveDestination_rois    = '/Users/jglendin/Desktop/HumanMSmasks/newROIs/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% define the file formats of images currently in 'input', and desired format of output
% INPUT ONLY can be a list of formats contained within a character array 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_rawdata_format            = {'.tiff','.tif'}; % can be a list

    segm_saveFMT_adjImage           = '.tiff'; 
    segm_saveFMT_foregroundMask     = '.png'; 

    reg_saveFMT_adjImage            = segm_saveFMT_adjImage;
    reg_saveFMT_foregroundMask      = segm_saveFMT_foregroundMask;

    chooseROI_saveFMT_allROIs       = '.tiff';
    chooseROIS_saveFMT_roiLocations = '.png'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do you want to process a particular image?
% Each script will default to looping through each and every image every the directory its
% been set to. You might want to override this... meaning you want to process only once.
% - you MUST change "loopThroughDirectory" to false! then indicate the sample ID and stain ID(s)
% - Ensure that there are not multiple files with the same sample ID, same stain
% ID, and matching extension (as indicated above as the input format). Multiple files
% matching these conditions will throw an error
% - For picking ROIs + segmentation, which operate w/1 file at a time, it will be assumed
% this file is described by stain1. Stain2 can have nothing or nonsense - it will be ignored until you try
% run a registration script.
% - For registration, assume that the "moving" image is stain2, and the "fixed" (ie the reference
% or template image) to which the moving image will be registered is stain1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doNOTloopThroughDirectory_justUseThis = false;

sampleID      = '4817'; % ==================================================== %
Stain1        = 'cd31'; % ==================================================== %
Stain2_moving = 'E7';   % ==================================================== %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Customizing function features and pararmeters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Segmentation --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a 1.25GB image takes ~1-2mins after it has been scaled down 2.8
        seg_scaleFactor = 2.8;

    % save scaled downimage? or try to resample image back to original size before save?
        seg_resizeBig = false;
    
    % after you have processed a raw file, would you like to delete it? If you are working
    % out of a folder of copies, it can help you stay organized...
        seg_deleteYES = true;
    
    % in the event that one file in your rawdata directory has a mask already in the save
    % folder, do you want to double check it or just skip?
        seg_doubleCheckSavedMasks = false;
    
    % optional - include a prefix identifier for saved files. It must start with an underscore!
         seg_versionID = '_v1';  %can also be left blank with char() 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Registration  --
    % Include optional step where a GUI opens to mediate manual coarse adjust rotation??
    % This is necessary if and only if your images are potentially > 90degrees unaligned
        reg_coarseAdj = true;

    % Include another optional step where again a GUI appears that allows you to, 
    % temporarily, manually adjust the shape of the binary mask of your foreground 
    % segmentation to make it more advantageous. My intended use case is if an image has a
    % unique rip, tear, debris, which will hinder its ability to 
        reg_manPtRepositioning = false;

    % scale factor #2. If your input to this function was scaled down for segmentation and
    % saved that way, then you can make this scale factor = 1
        reg_scaleFactor2 = 1.4; 
    
    % Set how many control points we assign to each image! In total we will have this
    % number*4 + 9 ...
        reg_numCPointsPerSide = 22;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choosing ROIs --
        % How many ROIs to pick?
        roi_numROIs = 10;

    % Size of ROI's? This is measured in pixels, and should be calculated based on the RAW
        % image. Give size in [Height pixels, Width pixels]:
        roi_sizeROI = [700, 700];

    % Compound masks? inlcude them here, and ROI will be sampled from their union. White should be foreground. 
        % TO OMIT - write either 0 or "none" below. 
        % TO INCLUDE masks, write the full path where additional masks will be. If multiple,
        % make a cell array
    roi_additionalMasks = 0; % or {'path/to/additionalmask1', 'path/to/additionalmask2',...}
% ==================================================== %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
% ==================================================== %
fileDirectories = struct('rawData', directOfImages,...
    'saveSegm_adjImages', saveSegm_adjImages,...
    'saveSegm_foregroundMask', saveSegm_foregroundMask, ...
    'saveRegistraion', saveRegistraion,...
    'saveDestination_rois', saveDestination_rois);

% ensure that each path in this array has a terminal '/'
fileDirectories = structfun(@(x) strcat(x, '/'), fileDirectories, 'UniformOutput', false);
fileDirectories = structfun(@(x) replace(x, '//', '/'), fileDirectories, 'UniformOutput', false);

fileFMTS = struct('rawDataFMT', [],...
    'segm_saveFMT_adjImage', segm_saveFMT_adjImage, ...
    'segm_saveFMT_foregroundMask', segm_saveFMT_foregroundMask,...
    'reg_saveFMT_adjImage', reg_saveFMT_adjImage,...
    'reg_saveFMT_foregroundMask', reg_saveFMT_foregroundMask,...
    'chooseROI_saveFMT_allROIs', chooseROI_saveFMT_allROIs,...
    'chooseROIS_saveFMT_roiLocations', chooseROIS_saveFMT_roiLocations);
fileFMTS(1).rawDataFMT = input_rawdata_format; %pulled out like this to be sure that cell array of formats doesn't fill in fileFMTS(2+3).rawdata
% ==================================================== %
pickME = struct('sampleID', [],'stain1', [], 'stain2', []);
if doNOTloopThroughDirectory_justUseThis
   pickME = struct('sampleID', sampleID,'stain1', Stain1, 'stain2', Stain2_moving);
end
% ==================================================== %
chooseROI = struct('numROIs', roi_numROIs, 'sizeROI', roi_sizeROI, 'roi_additionalMasks', roi_additionalMasks);
reg = struct('reg_coarseAdj',reg_coarseAdj, 'reg_manPtRepositioning', reg_manPtRepositioning,...
    'reg_scaleFactor2', reg_scaleFactor2, 'numCPointsPerSide', reg_numCPointsPerSide);
seg = struct('seg_scaleFactor', seg_scaleFactor, 'seg_resizeBig', seg_resizeBig, 'seg_deleteYES', seg_deleteYES, 'seg_doubleCheckSavedMasks', seg_doubleCheckSavedMasks);
savePrefs = struct('seg_versionID', seg_versionID);
% ==================================================== %
mySettings = struct('directories', fileDirectories,...
    'fileFormats', fileFMTS, 'seg', seg, 'chooseROI', chooseROI,...
    'reg', reg,'pickME', pickME, 'savePrefs', savePrefs);
% ==================================================== %


end

