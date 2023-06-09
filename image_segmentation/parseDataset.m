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

        %% OPT 1 -  we make the datastore containing specific images specified
        % if these fields are empty then we know that it was indicated in settings to process whole directory
        if ~isempty(mysettings_prefs.pickME.sampleID) && ~isempty(mysettings_prefs.pickME.stainID)
            sampleID = mysettings_prefs.pickME.sampleID;
            stainID = mysettings_prefs.pickME.stainID;

            % Loop over file extensions and add matches to Fnames
            Fnames = cell(1, numel(rawDataFormat));  % Initialize
            for i = 1:numel(rawDataFormat)
                rawFiles = dir(fullfile(rawDataLocation, ['*' rawDataFormat{i}]));
                Fnames{i} = {rawFiles.name};
            end
            Fnames = [Fnames{:}];

            %look for the files with names that match your Identifiers sample stain
            fIdx = contains(extractBefore(Fnames,'.'), sampleID,'IgnoreCase',true) & contains(extractBefore(Fnames,'.'), stainID,'IgnoreCase',true);

            if sum(fIdx)==0
                error(strcat('There is no file in this directory that shares the sampleID (', sampleID,'), stainID (', stainID,'), versionID and extension indicated'));
            elseif sum(fIdx)>1
                error(strcat('More than 1 image has been found within the raw data folder sharing the sampleID (', sampleID,'), stainID (', stainID,'), versionID, and extension you have indicated.'));
            end

            fileName = Fnames{fIdx};
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
        
        if class(data.imageNames) == 'char'
            nFiles = 1;
            data.imageNames = {data.imageNames};
        else
            nFiles = numel(data.imageNames);
        end

        for k = 1:nFiles

            % search for each part of filename within save_fL1 and save_fL2
            searchterms =  split(data.imageNames(k), '_');
            
            % fix the first term which will have the whole path
            extraSearchParts1 = split(searchterms{1}, '/');
            searchterms{1} = extraSearchParts1{numel(extraSearchParts1)};
            
            % fix the final term which will have the esxtension
            extraSearchParts2 = split(searchterms{3}, '.');
            searchterms{3} = extraSearchParts2{1};
            
            % This should just be sample ID, stain ID, and group ID, nothing else
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
            else
                error('nothing found at end of parser')
            end
        end

    case 'registration'

        %% Dropping soon


    case 'choosingROIs'
        if ~isempty(mysettings_prefs.pickME.sampleID) && ~isempty(mysettings_prefs.pickME.stain1)
            sampleID = mysettings_prefs.pickME.sampleID;
            stainID = mysettings_prefs.pickME.stain1;

            rgbIms = dir(fullfile(mysettings_prefs.directories.saveSegm_adjImages, ['*' mysettings_prefs.fileFormats.segm_saveFMT_adjImage]));
            Fnames1 = [rgbIms.name{:}];

            maskIms = dir(fullfile(mysettings_prefs.directories.saveSegm_foregroundMask, ['*' mysettings_prefs.fileFormats.segm_saveFMT_foregroundMask]));
            Fnames2 = [maskIms.name{:}];

            %look for the files with names that match your Identifiers sample stain
            fIdx = contains(Fnames1, sampleID,'IgnoreCase',true) & contains(Fnames1, stainID,'IgnoreCase',true);
            fIdx2 = contains(Fnames2, sampleID,'IgnoreCase',true) & contains(Fnames2, stainID,'IgnoreCase',true);

            if sum(fIdx)==0 || sum(fIdx2)==0
                error(strcat('You are missing at least 1 of the specified images. Specifications were: sampleID (', sampleID, '), stainID (', stainID, ')'));
            end
            if sum(fIdx)>1 || sum(fIdx2)>1
                % Regular expression pattern
                pattern = '_v(\d+).';
                versionNumbers1 = zeros(numel(rgbIms.name), 1, 'double');
                for k = 1:numel(rgbIms.name)
                    match1 = regexp(rgbIms.name(k), pattern, 'tokens', 'once');
                    if ~isempty(match1)
                        versionNumbers1(k, 1) = str2double(match1{1});
                    end
                end
                versionNumbers2 = zeros(numel(maskIms.name), 1, 'double');
                for b = 1:numel(maskIms.name)
                    match2 = regexp(maskIms.name(b), pattern, 'tokens', 'once');
                    if ~isempty(match2)
                        versionNumbers2(b, 1) = str2double(match2{1});
                    end
                end

            end
            [~, versionIdx1] = max(versionNumbers1);
            fileName1 = Fnames1{versionIdx1};
            [~, versionIdx2] = max(versionNumbers2);
            fileName2 = Fnames2{versionIdx2};

            data.imageDS = imageDatastore({fullfile(saveSegm_adjImages, fileName1), ...
                fullfile(saveSegm_foregroundMask, fileName2)});

            data.imageNames = data.imageDS.Files{1};


        else
            RGB_imageDatastore = imageDatastore(mysettings_prefs.directories.saveSegm_adjImages, 'FileExtensions', mysettings_prefs.fileFormats.segm_saveFMT_adjImage);
            fileList = RGB_imageDatastore.Files;
            [~, RGB_basefilenames, ~] = fileparts(fileList);

            nImages = numel(RGB_basefilenames);

            %before we can extract using underscores as landmarks, we must ensure that each file has
            %4! no more, no less. otherwise it will throw an error
            for h = 1:nImages
                currentFilename = RGB_basefilenames{h};
                underscores = strfind(currentFilename, '_');
                if numel(underscores)==5  % when this occurs defintiely we want to delete the 3rd
                    idx = underscores(3);
                    newFilename = strcat(currentFilename(1:idx-1), currentFilename(idx+1:end));
                    movefile(currentFilename, newFilename);
                    RGB_basefilenames{h} = newFilename;
                elseif numel(underscores)>5 || numel(underscores)<4
                    error('filenames are whack')
                end

            end

            % predefine this for latter
            pattern = '_v(\d+)';

            % from every file name, now we can extract the prefix (which contains the 3 key identifiers)
            allfileNameParts = split(RGB_basefilenames, '_');
            allRGB_prefix = strcat(allfileNameParts(:, 1),'_',allfileNameParts(:, 2), '_',allfileNameParts(:, 3));

            % It will happen that there are multiple files of the same image/,multiple image of same tissue
            uniquePrefixes = unique(allRGB_prefix);
            bestRGBpaths = cell(numel(uniquePrefixes),1); %loop through each unique prefix and decide on a best one
            for v = 1:numel(uniquePrefixes)
                prefix = uniquePrefixes(v);
                prefix_indices = contains(RGB_basefilenames, prefix);
                prefix_filenames = RGB_basefilenames(prefix_indices);
                version_numbers = regexp(prefix_filenames, pattern, 'tokens');
                version_numbers = cellfun(@(x) str2double(x{1}), version_numbers);
                [~, max_index] = maxk(version_numbers, 1);
                bestRGBpaths{v} = prefix_filenames{max_index};
            end
            unique_filepaths = fullfile(mysettings_prefs.directories.saveSegm_adjImages, strcat(bestRGBpaths,  mysettings_prefs.fileFormats.segm_saveFMT_adjImage));
            imds_uniqueRGB = imageDatastore(unique_filepaths);
            clear("RGB_imageDatastore");

            nImages = numel(imds_uniqueRGB.Files);
            maskFilenames{nImages} = '';
            for p = 1:nImages
                % start looping through the RGB images one at a time to find associated masks

                prefix = allRGB_prefix{p}; 

                % Search for the mask image filename that matches prefix
                maskIm_matches = dir(fullfile(mysettings_prefs.directories.saveSegm_foregroundMask, [prefix, '*' , mysettings_prefs.fileFormats.segm_saveFMT_foregroundMask]));

                if ~isempty(maskIm_matches)
                    if size(maskIm_matches.name, 1) == 1
                        maskFilenames{p} = fullfile(mysettings_prefs.directories.saveSegm_foregroundMask, maskIm_matches.name);
                    elseif numel(maskFilenames) > 1
                        matchesIMDS =imageDatastore(fullfile(mysettings_prefs.directories.saveSegm_foregroundMask, maskIm_matches.name));
                        matches = readall(matchesIMDS);
                        findVersionNum = regexp(matchesIMDS.Files(p), pattern, 'tokens', 'once');
                        maxVersionNum = find(findVersionNum== max(findVersionNum));
                        if numel(maxVersionNum)== 1
                            maskFilenames{p} = matchesIMDS.Files{maxVersionNum};
                        else
                            rgbIm = read(imds_uniqueRGB, p);
                            diceScore = zeros(numel(matches),1);
                            for r = 1:numel(matches)
                                diceScore(r,1) = dice(rgbIm, matches{r});
                            end
                            [~, bestIdx] = maxk(diceScore, 1);
                            maskFilenames{p} = matchesIMDS.Files{bestIdx};
                        end
                    else
                        disp('EEK');
                    end
                else
                    disp('EEK');
                end
            end

            data.maskIMDS = imageDatastore(maskFilenames);
            data.rgbIMDS = imds_uniqueRGB;

            [~, data.maskNames, ~] = fileparts(data.maskIMDS.Files);
            [~, data.rgbNames, ~] = fileparts(data.rgbIMDS.Files);


            % ====================================
        end

end

