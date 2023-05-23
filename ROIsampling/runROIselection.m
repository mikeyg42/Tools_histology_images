function runROIselection
% This is the principle  function one should call in order to initiate the selection of
% ROIs from a large, whole slide image. 

% Before rrunning this script, you will need to 1. modify the settings and preferences
% .m file entitled "setts_and_prefs.m", which will inform this code what it needs to do
% and where the data is to do so with. 2. you will also need to have already generated a
% binary iamge for each image you wish to select ROIs from, in which the foreground, must 
% be white, ie =1, atop a black background. 

% Note that even if you only want to run 1 image, this function will still work the same.
% The data parsing function will make a datastore of length 1 and eas soon as youre done
% it will save. 

% Do not forget to adjust your desired ROI size to be consistent with any preprocessing
% image resizing you may have included in your timeline. If you resized to 0.25x, then
% your 500x500 pixel ROI will "really" be selecting an area 2000x2000 pixel.

% The circle packing method implemented in the choosingROIs branch is but half completed for the time being!

% -Michael Glendinning, 2023

mySettings = setts_and_prefs;
data = parseDataset(mySettings, 'choosingROIs'); 



if contains(mySettings.chooseROI.roi_method, 'random', 'IgnoreCase',true)
    
    savePath = mySettings.directories.saveDestination_rois;
    if ~exist(fullfile(savePath, filesep, 'savedOutput'), 'file')
        mkdir(fullfile(savePath, filesep, 'savedOutput'));
    end
    savePath = fullfile(savePath, filesep, 'savedOutput');
    
    for nFile = 1:size(data.rgbIMDS.Files,1)
        %% run the script
        [all_roi_masks, array_of_individual_masks, array_of_ROIs] = userOnLoop_ROIselection(nFile);
        disp(strcat('done processing :',data.rgbNames(nFile, 1), ' beginning to save!'));

        %% save!!!
        % 1. Binary image that is the composite of all mask locations
        filename = strcat('All_ROImasks_', data.maskNames(nFile, 1), mySettings.fileFormats.chooseROIS_saveFMT_roiLocations);
        filepath_all = fullfile(savePath, filesep, filename);
        imwrite(all_roi_masks, char(filepath_all))

        % 2. Each individual ROI's Location AND 3. the actual cropped out ROIs
        filename_individ = strcat('individual_ROImasks_', data.maskNames(nFile, 1), mySettings.fileFormats.chooseROIS_saveFMT_roiLocations);
        filepath_individ = fullfile(savePath, filesep, filename_individ);
        filename_rois = strcat('individual_RGB_ROIs_', data.rgbNames(nFile, 1), mySettings.fileFormats.chooseROI_saveFMT_allROIs);
        filepath_rois = fullfile(savePath, filesep, filename_rois);
        for im  = 1:mySettings.chooseROI.numROIs
            if im == 1
                imwrite(array_of_individual_masks{im}, char(filepath_individ), 'WriteMode', 'overwrite');
                imwrite(array_of_ROIs{im}, char(filepath_rois), 'WriteMode', 'overwrite', 'Compression', 'none');
            else
                imwrite(array_of_individual_masks{im}, char(filepath_individ), 'WriteMode', 'append');
                imwrite(array_of_ROIs{im}, char(filepath_rois), 'WriteMode', 'append', 'Compression', 'none');
            end
        end
        disp(strcat('done saving :',data.rgbNames(nFile, 1) ));
    end

else
    disp(strcat('note: this script is incomplete still! Comment out the return in runROIsection.m',...
        'and try again if you wish to debug or help get me past a third of the way through'));
    return
    choosingROIs;
end

disp('Complete. Sucess.')