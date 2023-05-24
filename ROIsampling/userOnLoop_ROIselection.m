function [allROIsamples, FINmasks, extractedROIs] = userOnLoop_ROIselection(fileNumber)

mySettings = setts_and_prefs;

% call your parseDataset function, which interprets settings and loads up the dataset into
% a structure called "data"
data = parseDataset(mySettings, 'choosingROIs'); % inside data will be a image datastore from which you can read in images

% load up and check over your input files
imgCurrent = readimage(data.rgbIMDS, fileNumber);
imgCurrent = im2double(imgCurrent);
maskCurrent = logical(readimage(data.maskIMDS, fileNumber));
if size(imgCurrent, 1:2) ~= size(maskCurrent)
    maskCurrent = imresize(maskCurrent, size(imgCurrent,1:2), {@oscResampling, 4});
end

if ~strcmp(mySettings.chooseROI.roi_additionalMasks, 'none')
    for k = 1:numel(mySettings.chooseROI.roi_additionalMasks)
        extraMask = imread(mySettings.chooseROI.roi_additionalMasks{k});
        if size(extraMask) ~= size(imgCurrent, 1:2)
            extraMask = imresize(extraMask, size(imgCurrent,1:2), {@oscResampling, 4});
        end
        
        if ~islogical(extraMask)
            extraMask = logical(extraMask);
        end
        
        maskCurrent = maskCurrent & extraMask;
    end
end
close all force

% make sure imshow settings will give you largest image field possible 
iptsetpref("ImshowBorder","tight");

%% start of the processing to generate 10 ROI's
close all

[FINmasks, rectVerts, allROIsamples] = intializeRandomROIs(imgCurrent, maskCurrent, mySettings);

visualizeROIs(rectVerts, allROIsamples, imgCurrent, mySettings);

% open up choice of ROI replacement GUI   
CHOICE = showROIs_pickOne_GUI(mySettings);
while CHOICE<=mySettings.chooseROI.numROIs
% process user response:
        [allROIsamples,rectVerts, FINmasks] = replaceROI(CHOICE, rectVerts,FINmasks, allROIsamples, maskCurrent, imgCurrent, mySettings);  
        CHOICE = showROIs_pickOne_GUI(mySettings);
end
   
extractedROIs = extractROIs(imgCurrent, rectVerts, mySettings);

end

% =============== end principle code block =================================
% ==========================================================================
% =============== start initialization fcn =================================

function [FINmasks, rectVerts, allROIsamples] = intializeRandomROIs(imgCurrent, maskCurrent, mySettings)

%initialize 1. mask for storing distrution of ROIs map and 2. counter of ROIs 
% 3. cell arrays to store cartesian coordinates for ROIs, 4. each individial mask
allROIsamples= false(size(imgCurrent, 1:2));
numROIs = 1;
rectVerts = zeros(mySettings.chooseROI.numROIs, 5, 'double');
FINmasks{mySettings.chooseROI.numROIs} = [];

% pre calculate what you can:
tooManyROIs = (1+mySettings.chooseROI.numROIs);
areaROI= prod(mySettings.chooseROI.sizeROI);
imgSz = size(imgCurrent, 1:2);

while numROIs < tooManyROIs
%if mySettings.chooseROI.yesRotation

    [maskRectROI, rectVerts] = makeRandomWindowMask(imgSz, rectVerts, mySettings, numROIs);

    MASKintersect = maskRectROI & maskCurrent;
    overlappingROIs = maskRectROI & allROIsamples;

    if sum(MASKintersect(:)) == areaROI && sum(overlappingROIs(:)) == 0
        allROIsamples = allROIsamples | maskRectROI;
        FINmasks{numROIs} = maskRectROI;
        disp(strcat('done with number : ', num2str(numROIs)));
        
        numROIs = numROIs+1;
    end
end

end

% ================ end initialization fcn ==================================
% ==========================================================================
% ================ start random selection of ROI fcn =======================

function [maskRectROI, rectVerts] = makeRandomWindowMask(imgSz, rectVerts, mySettings, sampNum)
    % select a a random rectangular region using randomWindow2d builtin function
    win = randomWindow2d(imgSz(1:2), mySettings.chooseROI.sizeROI);
    newRectVerts = [win.XLimits(1), win.YLimits(1),diff(win.XLimits)+1, diff(win.YLimits)+1];
   
    % convert this limited Rectangle object, defined only by limits, into a noraml rectangle,
    % which we can importatnly specify to be rotatable if desired.
    rectROI = images.roi.Rectangle('Position', newRectVerts, 'Rotatable',true);
    
    if mySettings.chooseROI.roi_YesRotation
        rotAngle=randperm(360, 1);
        rectROI.RotationAngle = rotAngle;
        rectVerts(sampNum, 5) = rotAngle;
    end
    
    rectVerts(sampNum, 1:4) = rectROI.Position;
    maskRectROI = logical(poly2mask(rectROI.Vertices(:, 1), rectROI.Vertices(:,2), imgSz(1),imgSz(2)));

end 

% ================ end random selection of ROIs fcn ========================
% ==========================================================================
% ================ start visualization fnc =================================

function visualizeROIs(rectVerts, allROIsamples, imgCurrent, mySettings, varargin)

close all force

% Create a figure to visualize the image in. modify settings to ensure quality and
% efficiency are maximized
fHandle =figure;
set(fHandle,'doublebuffer','off'); %improves speed
ax1 = axes(fHandle);
imshow(imgCurrent, 'Parent', ax1);
set(ax1, 'xlimmode','manual',...
    'ylimmode','manual',...
    'zlimmode','manual',...
    'climmode','manual',...
    'alimmode','manual',...
    'Units','pixels');
fHandle.Visible = 'off';

% get the image extents from the axes properties
width = floor(diff(get(ax1,'XLim')));
height = floor(diff(get(ax1,'YLim')));

% parse whether or not a replacement was included in the file calll or not.
if all(size(varargin, 1:2)>0)
    if ndims(varargin{1}) == 2
    replacementMask_aData = logical(imresize(varargin{1} , [height, width], {@OSCresampling, 4}));
    
    assert(isreal(varargin{2}))
    numberMagenta = varargin{2}; 

    else % ie if its nonempty but not a gray image....
        error('visualizeROI function: varargin{1} is a unusable format we cannot handle. kindly fix it before re-trying to run')
    end
else
    replacementMask_aData = false([height, width]); 
    if any(isnan(rectVerts))
        error('visualizeROI function expected that in the event that there are rectVerts=NaN, that varargin{1} and {2} would contain information of to replacement those NaNs...')
    end
end

nROI = mySettings.chooseROI.numROIs;
sizeROI = mySettings.chooseROI.sizeROI;

% convert rectangle coordinates to centoids, adding some displacement to facilitate visualization:
halfxy = floor(rectVerts(:, 3:4)./2);
offset4 = round(sizeROI(1:2)/4);
displacement = repmat(offset4, size(halfxy,1), 1);
textLocation = rectVerts(:,1:2)+halfxy-displacement;

% prepare flat blocks of opaque colors which will become alpha data showing ROI in situ
% locatio 
ONES = mat2gray(ones([height, width], 'double'));
ZEROS = mat2gray(zeros([height, width], 'double'));
compositeOverlay = cat(3, ZEROS, ONES, ZEROS);
magentaOverlay = cat(3, ONES, ZEROS, ONES);
% merge the images using replacementMask_aData as a guide
compositeOverlay(repmat(replacementMask_aData, [1 1 3])) = magentaOverlay(repmat(replacementMask_aData, [1 1 3]));

hold on % now show overlay the greeen (maybe magenta too) brightly color block
gr_mg_Overlay = imshow(compositeOverlay);
hold off

% Use the allSamples image (merged with the replacement ROI mask) to be the alphaDataMap
alphaData_mask = imresize(allROIsamples, [height width], {@oscResampling, 4});
alphaData_mask = alphaData_mask | replacementMask_aData;

% Set the alphaDataMask so that green/magenta only is visible above the ROIs selected
alphaData_mask = im2double(alphaData_mask | replacementMask_aData);
alphaData_mask = alphaData_mask.*0.4; 
set(gr_mg_Overlay, 'AlphaData', alphaData_mask)

% Label these ROIs with their numeric identification
for jj = 1:nROI        
        text(textLocation(jj,1), textLocation(jj,2), num2str(jj), 'FontSize', 14, 'FontWeight', 'Bold');
end

drawnow;
fHandle.Visible = 'on';
end

% ================= end visualization fcn ==================================
% ==========================================================================
% ================= start replacing ROI fcn ================================

function [allROIsamples,rectVerts, FINmasks] = replaceROI(CHOICE, rectVerts,FINmasks, allROIsamples, maskCurrent, imgCurrent, mySettings)

close all

allROIsamples(FINmasks{CHOICE}) = false;
FINmasks{CHOICE} = false(size(allROIsamples));
rectVerts(CHOICE, 1:5)= [NaN, NaN, NaN, NaN, NaN];

sizeROI = mySettings.chooseROI.sizeROI;
imgSz = size(maskCurrent);

while true

    [newROImask, newRectVerts] = makeRandomWindowMask(imgSz, rectVerts, mySettings, CHOICE);

    crossOutsideBounds_test= newROImask & maskCurrent;
    overlapAnotherROI_test = newROImask & allROIsamples;
    
    if  sum(newROImask(:)) == 0 || sum(crossOutsideBounds_test(:)) ~= sum(newROImask(:)) || sum(overlapAnotherROI_test(:)) > 0
        continue
    else

        rectVerts(CHOICE, 1:5) = newRectVerts(CHOICE, 1:5);

        visualizeROIs(rectVerts, allROIsamples, imgCurrent, mySettings, newROImask, CHOICE);
        
        selection = questdlg('Evaluate new ROI', 'ROI inquiry',...
            'ROI is good','ROI is bad, plz remake it', 'ROI is bad, plz remake it');

        if contains(selection, 'good')
            allROIsamples = allROIsamples | newROImask;
            FINmasks{CHOICE} = FINmasks{CHOICE} | newROImask;
            break

        end %end of if contains good block
    end %if block

end % end of the while loop

end

% ================= end replacing ROI fcn ==================================
% ==========================================================================
% ================= start selection GUI fcn ================================

function CHOICE = showROIs_pickOne_GUI(mySettings)
%% we ask the User here if they want to exit and move to another image, or replace an image that 
%% captures a a suboptimal tissue area (e.g. a rip, a hole)?

% We make a pop-up GUI that will prompt user with a list of options, which takes the form of
% a numbered list (the numbers on the list cooresponding to the in situ numbered tiles overlaid in
% the visualization that will also live on screen. The GUI will also have an option to
% say, "okay, please exit".

% Prepare the menu options for the GUI
numlist = num2cell(1:mySettings.chooseROI.numROIs);
exitChoice = {'all good!'};

% Instantiate the choice menu GUI 
CHOICE = menu('Replace any of these any ROIs?', [numlist, exitChoice]);
end

% ================= end selection GUI fcn ==================================
% ==========================================================================
% ================= start extract ROIs fcn =================================
function extractedROIs = extractROIs(imgCurrent, rectVerts, mySettings)


% Initialize variables
numROIs = size(rectVerts, 1);
extractedROIs = cell(numROIs, 1);

for im = 1:numROIs
    % Extract coordinates and rotation angle for the current ROI
    rectCoords = rectVerts(im, 1:4);
    rotationAngle = rectVerts(im, 5);
    
    % Calculate center point of the unrotated image
    centerX = rectCoords(1) + rectCoords(3)/2;
    centerY = rectCoords(2) + rectCoords(4)/2;
    translation1 = [1 0 -centerX; 0 1 -centerY; 0 0 1];

    % Rotate the cropped image
    rotation = [cosd(rotationAngle) -sind(rotationAngle) 0; sind(rotationAngle) cosd(rotationAngle) 0; 0 0 1];

    % invert the initial translation
    inverseTranslation = [1 0 centerX; 0 1 centerY; 0 0 1];

    % combine these all together and apply it to the image:
     transformationMatrix = inverseTranslation * rotation * translation1;
     tform = affinetform2d(transformationMatrix);
    
    % apply the transformation insisting the ouitput view is unchanged from the input
    sameAsInput = affineOutputView(size(imgCurrent),tform, 'BoundsStyle','SameAsInput');
    transformedImage = imwarp(imgCurrent, tform, 'cubic', 'OutputView', sameAsInput);

    % Crop the image and store in the cell array you preallocated 
    rect = images.spatialref.Rectangle([rectCoords(1), rectCoords(1)+rectCoords(3)-1],...
       [rectCoords(2), rectCoords(2)+rectCoords(4)-1]); % subtract to accomodate 1-pixel indexing, I think? 

     extractedROIs{im} = rot90(imcrop(transformedImage, rect), -1); 
end


end

