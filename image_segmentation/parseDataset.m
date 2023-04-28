function data = parseDataset(mysettings_prefs, SCRIPT)
%mysettings_prefs = struct('directories', fileDirectories, ...
%    'fileFormats', fileFMTS, 'seg', seg, 'chooseROI', chooseROI,...
%    'reg', reg,'pickME', pickME);
% Michael Glendinning, 2023

switch SCRIPT
    case 'segmentation'
        % preallocate structural array called data
        data = struct('imageDS', [], 'imageNames', [], 'segInfo', []);
        
        % two possibilities, running every file in the directory or picking 1 specific file
        rawDataFormat = mysettings_prefs.fileFormats(1).rawDataFMT;
        rawDataLocation = mysettings_prefs.directories.rawData;
        
        %% OPT 1 -  we make the datastore containing every image in directory
        % if these fields are empty then we know that it was indicated in settings to process whole directory
        if ~isempty(mysettings_prefs.pickME.sampleID) && ~isempty(mysettings_prefs.pickME.stain1)
            sampleID = mysettings_prefs.pickME.sampleID;
            stainID = mysettings_prefs.pickME.stain1;
            
            % Loop over file extensions and add matches to Fnames
            Fnames = cell(1, numel(rawDataFormat));  % Initialize
            for i = 1:numel(rawDataFormat)
                rawFiles = dir(fullfile(rawDataLocation, ['*' rawDataFormat{i}]));
                Fnames{i} = {rawFiles.name};
            end
            Fnames = [Fnames{:}];
            
            %look for the files with names that match your Identifiers sample stain
            fIdx = contains(Fnames, sampleID,'IgnoreCase',true) & contains(Fnames, stainID,'IgnoreCase',true);
            
            if sum(fIdx)==0
                error(strcat('There is no file in this directory that shares the sampleID (', sampleID,'), stainID (', stainID,'), versionID and extension indicated'));
            elseif sum(fIdx)>1
                error(strcat('More than 1 image has been found within the raw data folder sharing the sampleID (', sampleID,'), stainID (', stainID,'), versionID, and extension you have indicated.'));
            end
            
            fileName = Fnames{find(fIdx)};
            data.imageDS = imageDatastore(fullfile(rawDataLocation, fileName));
            data.imageNames = data.imageDS.Files{1};
        else
            %% OPT 2 -  we make the datastore containing every image in directory
            my_imageDatastore = imageDatastore(rawDataLocation, 'FileExtensions', rawDataFormat);
            fileList = my_imageDatastore.Files;
            fileList = extractBefore(extractAfter(fileList, rawDataLocation), '.'); %get only the final bit of filename after the parent directory, sans any file extension
            
            data.imageNames = fileList;
            data.imageDS = my_imageDatastore;
        end
        
        % CHECK SAVED DIRS for already finished work..Get list of saved adjRGB images and a list of saved Masks
        save_fL1 = string({dir(fullfile(mysettings_prefs.directories.saveSegm_adjImages, '*')).name});
        save_fL2 = string({dir(fullfile(mysettings_prefs.directories.saveSegm_foregroundMask, '*')).name});
        
        if ~isempty(mysettings_prefs.savePrefs.seg_versionID)
            vID = true;
            % get rid of any underscores
            versionID = setdiff(mysettings_prefs.savePrefs.seg_versionID, '_', 'stable');
        end
        
        for k = 1:numel(fileList)
            
            % search for each part of filename within save_fL1 and save_fL2
            searchterms =  split(data.imageNames(k), '_');
            myquery = {searchterms{1}, searchterms{2}, searchterms{3}};
            if ~vID
                matches1 = cellfun(@(x) contains(save_fL1', x,'IgnoreCase',true), myquery , 'UniformOutput', false);
                matches2 = cellfun(@(x) contains(save_fL2', x,'IgnoreCase',true), myquery, 'UniformOutput', false);
            else
                matches1 = cellfun(@(x) contains(save_fL1', x,'IgnoreCase',true), [myquery, {versionID}], 'UniformOutput', false);
                matches2 = cellfun(@(x) contains(save_fL2', x,'IgnoreCase',true), [myquery, {versionID}], 'UniformOutput', false);
            end
 % matches1/2 will both contain 3 or 4 vectors w/ same length as save_filelist1/2 
 % which will be all filled with logical values indicated string pattern matches. 
 % patterns are the 3 or 4 file naming convention componentswe searched, saved as 'myquery'
 
 % Cell2Mat results in 3-4 vectors being horzcat. Then, Prod in the 2nd dim propogates
 % an @AND across each row, s.t. a 1 in combined array is a complete match to each aspect
 % of myquery
 
            combinedArray1 = logical(prod(cell2mat(matches1), 2)); % adjusted RGB images
            combinedArray2 = logical(prod(cell2mat(matches2), 2)); % MASKS
            
            if sum(combinedArray1 | combinedArray2) == 0
                % pass this informatino along to the calling function so it knows what processing needs to be done
                data.segInfo(k).needAdjRGB = true;
                data.segInfo(k).needMask = true;
                
            elseif sum(combinedArray1) > 0 && sum(combinedArray2) == 0
                data.segInfo(k).needAdjRGB = false;
                data.segInfo(k).needMask =true;
                
            elseif sum(combinedArray1) == 0 && sum(combinedArray2) > 0
                data.segInfo(k).needAdjRGB = true;
                data.segInfo(k).needMask = false;
                
            elseif sum(combinedArray1) > 0 && sum(combinedArray2) > 0
                data.segInfo(k).needAdjRGB = false;
                data.segInfo(k).needMask = false;
            end
        end
        
    case 'registration'
        
    case 'chooseROIs'


% ====================================
end

end

