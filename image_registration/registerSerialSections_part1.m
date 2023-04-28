function registerSerialSections_part1
% This script takes as input 2 images along with a binary "masks" (where white indicates
% foreground and black indicated background) for each, and will attempt to register them
% semi-automatically using affine transformations. At the end of the script, part 2 is
% automatically initiated, which further improves the registration with first a local
% genometric transformation and then a diffeomorphic demons algorithm.


% Note: The script, as is, relies heavily upon your file naming convention.
%   Prepending my convention with your will prevent any issues. 
%   my convention requires three, underscore_seperated identifiers: 

%               {1} _ {2} _ {3} _ {optional 4}.tiff
%        ie. "GROUPid_SAMPLEid_STAINid_ExtraInfo.tiff   

% {1} - experimental group ID. I often use: "SICK" or "HEALTHY". 
%       
% {2} - SAMPLE id. I often use A, B, C, ...

% {3} - identification of the particular staining performed in that image
%       (H+E, ect.) For IHC usually the particular target antigen for
%       IHC-Dab (e.g. CD31, Sox17, ...)

% {4} - whatever other notes/identifiers you want. Your naming convention 
%       can just be here with "{1}_{2}_{3}_" tacked on front, and it will not 
%       mess up anything. Two exceptions:
%           - no additional underscores. there must be 3 and only 3 in the filename
%       	- do not repeat identical text thats in {1}, {2}, or {3}. 
%      e.g. [DO NOT EMULATE THIS]: '.../control_lumbar02_CD68_mybestCD68ControlSC.tiff'  <- BAD


% change these variable assignments! 
%  ||    ||
%  vv    vv     directorytosaveinto should be where you have your 
% ======================================================================================== %
directorywithimagefiles = '/Users/jglendin/Desktop/tiled_images_as_tiffs/';  %
destinationToSaveInto = '/Users/jglendin/Dropbox/Dropbox - Michael/processedMSimages';%
ID = 'A'; %===========================================================================  % 
Stain1_fixed = 'CD31'; %=================================================================  %
Stain2_moving = 'MCAM'; %=================================================================  %
% =======================================================================================  %
               
                %% ...THE GAME IS AFOOT %%
    
%% Part 0. Preprocessing steps 


%% preprocessing step 1. Load the image files 
% We need to load into the workspace the two images we will register. We
% also need masks delineating background/foreground. 
imDS = createDatastoreWithAll4Images(directorywithimagefiles, ID, Stain1_fixed, Stain2_moving);
myimages = cellfun(@(X)imresize(X, 0.5,{@oscResampling, 4}), readall(imDS), 'UniformOutput', false);

% ref_IMG will be also called "FIXED" to remind us to not perform any spatial transformations
ref_IMGadj = ensureDoubleScaled(myimages{2},true);
ref_IMGmask = logical(myimages{1});

% MOVING also be called "MOVING" because that is what we will
% morpho, maintaining the integrity of the stain we will later quantify
MOVING_pre = ensureDoubleScaled(myimages{4},true);
MOVINGmask_pre = logical(myimages{3});
clear myimages


%% preprocessing step 2. optional re-orientation GUI for coarse, manual oriention
close all force
pause(1);

% provide a visual on which to base decision of option GUI 
movingOverlay = imfuse(rgb2gray(MOVING_pre), MOVINGmask_pre, 'falsecolor', 'Scaling', 'none', 'ColorChannels', [2,1,2]);
fixedOverlay = imfuse(rgb2gray(ref_IMGadj), ref_IMGmask, 'falsecolor', 'Scaling', 'none', 'ColorChannels', [1,2,2]);
fOverlay = figure('Visible', 'on', 'HitTest', 'off'); 
axOverlay = axes(fOverlay, 'TitleFontSizeMultiplier', 1.4);
drawnow expose 
axOverlay.Title.String = 'Left = moving, Right = Fixed';
imshowpair(movingOverlay, fixedOverlay, 'montage', 'Parent', axOverlay);

answer1 = questdlg('Would you say that the moving and fixed images are roughly oriented the same??',...
    'coarse adjust orientation?','Yes', 'No, we must fix this', 'No, we must fix this'); 

% parse response
if contains(answer1, 'No') % this below calls the GUI
    [theta_cw, flip_LR] = RGB_reorientationGUI(MOVING_pre,ref_IMGadj);
    NquarterTurnsCW = theta_cw/90;
    fOverlay.Visible = 'off';
    % GUi gives you how many rotations to do... so do them!
    if flip_LR == 1 %if mirror image is req'ed
        MOVING_pre = fliplr(MOVING_pre);
        MOVINGmask_pre = fliplr(MOVINGmask_pre);
    end
    
    MOVINGadj = rot90(MOVING_pre, NquarterTurnsCW);
    MOVINGmask = rot90(MOVINGmask_pre, NquarterTurnsCW);
else
    fOverlay.Visible = 'off';
    MOVINGadj = MOVING_pre;
    MOVINGmask = MOVINGmask_pre;
end


%% preprocessing step 3. image enhancing
% backgrund correction and convert to grayscale WEIGHTS OF GRADIENTS -
% this is inversly related to gradient, s.t. smooth regions (tiny gradient)=big weights (i.e. white)
% and large gradients such as edges have small weight (black color).

% Background correction: (using mask and an extensively blurred and brightend image background subtraction).
se = strel('disk',25);

tophatFixed = 1 - imtophat(imcomplement(ref_IMGadj), se);
tophatMoving = 1 - imtophat(imcomplement(MOVINGadj), se);

% convert images into their WEIGHTED GRADIENTS :
% Weighting gradients are generally very FLAT. The algorithm
% prioritizes important original pixels like strong edges are kept
% unchanged. (its a 2nd deriv ultimately)

ref_IMGgray = rgb2gray(gradientweight(tophatFixed, 1.8, 'RolloffFactor', 1.25, 'WeightCutoff', 0.01));
MOVINGgray = rgb2gray(gradientweight(tophatMoving, 1.8, 'RolloffFactor', 1.25, 'WeightCutoff', 0.01));

%visualize
fOverlay.Visible = 'on';

movingOverlay = imfuse(MOVINGgray, MOVINGmask, 'falsecolor', 'Scaling', 'none', 'ColorChannels', [2,1,2]);
fixedOverlay = imfuse(ref_IMGgray, ref_IMGmask, 'falsecolor', 'Scaling', 'none', 'ColorChannels', [1,2,2]);
imshowpair(movingOverlay, fixedOverlay, 'montage', 'Parent', axOverlay);

%% step 4. optional GUI to optionally adjust MASK manually.
% ask about engaging in the option extra step to visualize the vertexes and
% hand arrange them to a certain extent.

answeringBox = questdlg('Would you like to modify either of the masks manually to facilitate their registration? A great option for handling small rips or folds',...
    'Tweak Any Vertecies manually?','Yes, please', 'No', 'Yes, please'); 

close 1

% parse response
if contains(answeringBox, 'Yes')
    MOVINGmask  = vertexTweaking_handdrawn(MOVINGmask,MOVINGgray, 120);
    ref_IMGmask  = vertexTweaking_handdrawn(ref_IMGmask,ref_IMGgray, 120);
end

%% step 5. background filter! do this by blurring background lots but not touching foreground using mask

% 0th. temporarily shrink
% 1st. blurring nonlinear Median Filter. 2nd. blurring with a gaussian filter, and 
% 3rd. H-minima transform is taken for each channel of the RGB image seperately. 
% 4th. bring back to the size before 0th step
% ~~ first the fixed image
tophatFixed_small = imresize(im2double(tophatFixed), round(size(tophatFixed, 1:2).*0.17), 'bilinear');
blurry1 = imgaussfilt(mikeMedianFilter(tophatFixed_small, 3, 1250, 'RGB'), 13, 'FilterSize', 25, 'FilterDomain', 'frequency'); 
blurry1 = cat(3, imhmin(blurry1(:,:,1),0.2), imhmin(blurry1(:,:,2),0.2), imhmin(blurry1(:,:,3),0.2)); %supress minima below 0.2
blurryFixed = imresize(blurry1, size(tophatFixed, 1:2), 'bilinear');
% ~~ repeat with the moving
tophatMoving_small = imresize(im2double(tophatMoving), round(size(tophatMoving, 1:2).*0.17), 'bilinear');
blurry2 = imgaussfilt(mikeMedianFilter(tophatMoving_small, 3, 1250, 'RGB'), 13, 'FilterSize', 25, 'FilterDomain', 'frequency'); 
blurry2 = cat(3, imhmin(blurry2(:,:,1),0.2), imhmin(blurry2(:,:,2),0.2), imhmin(blurry2(:,:,3),0.2)); %supress minima below 0.2
blurryMoving = imresize(blurry2, size(tophatMoving, 1:2), 'bilinear');

% ONLY use that super super blurry image on your BACKGROUND! Keep foreground pristine
ref_IMGgray(~ref_IMGmask) = blurryFixed(~ref_IMGmask);
MOVINGgray(~MOVINGmask) = blurryMoving(~MOVINGmask);

%% preprocessing step 6: image padding to preserve edges and prevent clipping from rotation
% also need to provide enough space for clipping that might occur after scaling and rotating!

padAmt = 200; %if substantial rotation will be required (>30deg), or images > 1.5 gigabyte add more pad.

[MOVINGmask, MOVINGgray, ref_IMGmask, ref_IMGgray] = pad4Images(padAmt, MOVINGmask, MOVINGgray, ref_IMGmask, ref_IMGgray);

%% DONE pre-processing

% recall the plan: 
% part 1. is rigid transformation (first rotation/scaling, then translation. 
% part 2. we set control points, then do nonlinear geometric transformation. 
% part 3. we do a diffeomorphic demons transformation. 

%% START ~ Part 1! %Rigid Geometric Transformations

%% Part 1. step 1: select 4 "corner points" around the sections' edge which will define rotation

% OPEN GUI, make sure we have each corner precisely right!
[movingMat, fixedMat, ~, ~] = selectCornersGUI(MOVINGmask, ref_IMGmask, MOVINGgray,ref_IMGgray);
% both of these output matrices have 4 rows cooresponding to 2 points. They
% are saved as [

%% ~~~~~~~~~~~~Rotation! and Scale!~~~~~~~~~~~~~~~~~~

% ============= estimate k 
k = estimateScaling(fixedMat, movingMat);

% ============= RESIZE HERE
MOVINGmask1 = imresize(MOVINGmask, k, {@oscResampling, 4});
MOVINGgray1 = imresize(MOVINGgray, k, {@oscResampling, 4});
movingMat2prime = movingMat.*k;

% ============= estimate rotation
theta = estimateRotation(movingMat2prime,fixedMat);

% ============= ROTATE HERE
[MOVINGgray2, MOVINGmask2, MOVINGrotatedMat] = execute3Rotation(MOVINGgray1, MOVINGmask1, movingMat2prime, theta);


% =%==%==%==%==%== a quick =%==%==%==%==%==%==%==%==%==
% ==%==%==%==% centroid detour! ==%==%==%==%==%==%==%==
% I want to include the centroid of each blob before next step. 
% It will help guide the translation substantially. Not only is it 125%
% more alignment data, to ensure the images are actually "on top of
% one another", but with very misaligned corners the forms can "veer", but
% the centroid grounds them by nature of its general invariance. 

ss = bwconncomp(MOVINGmask2, 8);
movingTable = regionprops('table', ss, 'Centroid', 'Area'); 
if size(movingTable, 2) > 1
toDelete = movingTable.Area<max(movingTable.Area);
movingTable(toDelete, :) = [];
elseif   size(movingTable, 2) < 1
    disp('ERROR');
end

MOVING_Mat(1:5, 1:2) = [MOVINGrotatedMat; movingTable.Centroid];
    
tt = bwconncomp(ref_IMGmask, 8);
fixedTable = regionprops('table', tt, 'Centroid', 'Area');
if size(fixedTable, 2) > 1
toDelete2 = fixedTable.Area<max(fixedTable.Area);
fixedTable(toDelete2, :) = [];
elseif  size(fixedTable, 2) < 1
    disp('ERROR');
end

fixedMat(1:5, 1:2) = [fixedMat; fixedTable.Centroid];   

% the MOVING mask has some numerical precision issues during the rotation
% that this corrects:
MOVINGmask2(MOVINGmask2>0.5) = 1; ref_IMGmask(ref_IMGmask>0.5) = 1;
MOVINGmask2(MOVINGmask2<0.5) = 0; ref_IMGmask(ref_IMGmask<0.5) = 0;

% ==%==%==%==% END Centroid Detour! =%==%==%==%= 
% %==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%=


%% ============= TRANSLATION! =============
% this fcn will also resize our images so that they are all the exact same square size
myOLDimages = {MOVINGmask2, MOVINGgray2, ref_IMGmask, ref_IMGgray}; 
[fixedMat_new,movingMat_new, myNEWimages] = preReg_padding_center(fixedMat, MOVING_Mat, myOLDimages);

% I padded with the wrong color somewhere - this is a lazy fix
im2 = myNEWimages{2}; myNEWimages{2} = imcomplement(imclearborder(imcomplement(im2)));
im4 = myNEWimages{4}; myNEWimages{4} = imcomplement(imclearborder(imcomplement(im4)));

%figure; montage(myNEWimages);

%% save your work here
part1 = struct([]);
part1(1).fixedMat = fixedMat_new;
part1(1).movingMat = movingMat_new;
part1(1).images_part1 = myNEWimages;
save(handy.savePath, '-struct', 'part1', '-v7.3', '-nocompression');

%% proceed to ~part 2~
% ==%==%==%==% Leave this Script and go to PART 2 =%==%==%==%==%==%===%==%==

[D, tform, m_im, movedMask] = reigsterSerialSections_part2_nonrigid(fixedMat_new,movingMat_new, myNEWimages);

% %==%==%==%==%= Welcome Back! and Good Job! =%==%==%==%==%==%==%==%==%==%==%=
%% append this second data structure to the .mat file where you saved part 1 info

part2 = struct([]);

part2(1).displacementField = D;
part2(1).transformation = tform;
part2(1).movedIm = m_im;
part2(1).movedMask = movedMask;
save(handy.savePath, 'data2', '-v7.3', '-nocompression', '-append');

end

