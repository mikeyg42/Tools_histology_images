function [D, tform, movedIm_partway, movedMask] = registerSerialSections_part2_nonrigid(varargin)
%syntax: [D, tform, M_Im, movedMask] = registerSerialSections_part2_nonrigid(varargin)

%% Part 2 of registration: nonrigid geometric transformation with control points placed EXCLUSIVELY programically .
% name of the game is placing control points programmically, and with great care.

%Using these points, I try 3 different nonrigid registration (LWH /2deg+3deg polynomial) techniques.
% After user picks the best of those, I follow up with a quick call of Thirion's demons
% algorithm, after which, images are usually very much aligned smoothly. registration is
% applied simultaneously to the image and its mask.

%Michael Glendinning, 2023

tic;

id = 'MATLAB:polyshape:repairedBySimplify';
warning('off', id);

%% load up variables saved
if nargin ~= 3
   ld = load('/Users/jglendin/Dropbox - Michael/Dropbox/human ms code/my code/registration_data_B_EGFL7_PLP.mat',...
        '-mat');
    % ld = load('/your/Dir/registration_data_sampleID_stain1_stain2.mat',...
     %   '-mat');
    myNEWimages = ld.images_part1;
    fixedMat = ld.fixedMat;
    movingMat = ld.movingMat;
else
    myNEWimages = varargin{3};
    fixedMat = varargin{1};
    movingMat = varargin{2};
end

MOVING_gray = myNEWimages{2};
MOVING_mask = myNEWimages{1};
IMG_gray = myNEWimages{4};
IMG_mask = myNEWimages{3};

MOVING_gray(~MOVING_mask) = 1;
IMG_gray(~IMG_mask) = 1;

%% define evenly spaced points along each of the four edges -- will become the majority of control points

%----------------------------------------------------------------------------
% start first with the fixed image
%----------------------------------------------------------------------------

B = bwboundaries(IMG_mask, 'noholes');
bxp = B{1}(:,2);
byp = B{1}(:,1);

pgon_fixed_nosimp = polyshape(bxp, byp, 'SolidBoundaryOrientation', 'cw'); % "cw" for clockwise
pgon_fixed = rmholes(pgon_fixed_nosimp);
%pgon_fixed = simplify(pgon_fixed_nosimp,'KeepCollinearPoints',1);
%----------------------------------------------------------------------------
% ...now do the moving image
%----------------------------------------------------------------------------

A = bwboundaries(MOVING_mask, 'noholes');
axp = A{1}(:,2);
ayp = A{1}(:,1);

pgon_move_nosimp = polyshape(axp, ayp, 'SolidBoundaryOrientation', 'cw');
pgon_move = rmholes(pgon_move_nosimp);
%pgon_move = simplify(pgon_move_nosimp,'KeepCollinearPoints',1);


%-----------
% QUERY your turning distance to evaluate if everything is working...
% td = turningdist(pgon_move, pgon_fixed);
% if td <0.4
%     disp('polygons are VERY different, you might have an issue?');
% end

%----------------------------------------------------------------------------
%% We need to adjust the corner points because they need to be perfect and the past few steps might've caused minor shifts

%movingMat/FixedMat has the cornerpoints, and pgon is the most recent pgon
cornerIndxMoving = nearestvertex(pgon_move, movingMat(1:4,:));
cornerIndxFixed = nearestvertex(pgon_fixed, fixedMat(1:4,:));

%----------------------------------------------------------------------------
% we are now subdividing the edges of each poygon into 4 contiguous pieces
% (delineaated by our favorite corner points, naturally).
%----------------------------------------------------------------------------

% the astute observer might notice that the cumulative length of these arrays of vertices
% is four greater than the initial array. this is because each array will have one of our
% corners doubled up. we will rectify this later on, but for now its a good thing.

coordXY_pgon_move_splitsides(1).points = [pgon_move.Vertices(cornerIndxMoving(2, :):end,:);pgon_move.Vertices(1: cornerIndxMoving(1, :),:)];
coordXY_pgon_move_splitsides(2).points = pgon_move.Vertices(cornerIndxMoving(1, :): cornerIndxMoving(4, :),:);
coordXY_pgon_move_splitsides(3).points = pgon_move.Vertices(cornerIndxMoving(4, :): cornerIndxMoving(3, :),:);
coordXY_pgon_move_splitsides(4).points = pgon_move.Vertices(cornerIndxMoving(3, :): cornerIndxMoving(2, :),:);

coordXY_pgon_fixed_splitsides(1).points = [pgon_fixed.Vertices(cornerIndxFixed(2, :):end,:);pgon_fixed.Vertices(1: cornerIndxFixed(1, :),:)];
coordXY_pgon_fixed_splitsides(2).points = pgon_fixed.Vertices(cornerIndxFixed(1, :): cornerIndxFixed(4, :),:);
coordXY_pgon_fixed_splitsides(3).points = pgon_fixed.Vertices(cornerIndxFixed(4, :): cornerIndxFixed(3, :),:);
coordXY_pgon_fixed_splitsides(4).points = pgon_fixed.Vertices(cornerIndxFixed(3, :): cornerIndxFixed(2, :),:);

%% Use curve fitting to place a LOT of perfect, amazing, control points
%% (i.e. control points programmically positioned evenly along the entire edge of the tissue.
% Because we've already coarse aligned the images, these points should be well
% aligned across images...


%------- SET NUM POINTS ALONG AN EDGE---------------------------------------------------------------------
pointsPerSide = 20;

% this value x4 and +9 (to account for 3x3 grid of middle points),
% this sum will equal our maximum # control points. (e.g. if 15 points per side, we will
% have 69 total control points)
%----------------------------------------------------------------------------

%%
%preallocate for the big loop
evenlydistMoving = zeros(pointsPerSide*4,2, 'double');
evenlydistFixed = zeros(pointsPerSide*4,2,'double');
midGrid_m{4} = zeros(3,2);
midGrid_f{4} = zeros(3,2);

%% ||~~~~-~~~~||~~~~-~~~~|| START of LOOP ||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||
%%      ||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||

counter = 1;
for cornerN = 1:4
    
    PointsToFit_move =  coordXY_pgon_move_splitsides(cornerN).points;
    PointsToFit_fix = coordXY_pgon_fixed_splitsides(cornerN).points;
    
    %~~~~-~~~~||~~~~-~~~~||~call the curve fitting script!||~~~~-~~~~||~~~~-~~~~||
    
[pointData]  = curveFittingOfTissueBorders(pointsPerSide,cornerN, PointsToFit_move, PointsToFit_fix, MOVING_gray, IMG_gray);

    % Q: reorder the midgrid points? ~~~~~~~-~~~~~~~~
    %  NOPE! in old versions of the code I had to flip some of these coordinates. But in
    %  its current form the points should be already arranged (top of matrix to bottom like so:
    %         1  2  3 
    %         |  |  |
    %     3 --+--+--+-- 3
    %     2 --+--+--+-- 2
    %     1 --+--+--+-- 1
    %         |  |  |
    %         1  2  3
    %~~~~-~~~~~~~~-~~~~~~~~-~~~~~~~~--~~~~~~~~-~~~~~~~~-~~~~~~~~-~~~~~~~~~~-~~~~~~~~~~-
    % parsing curve fitting script results:
    %the results of that script are a structural array (admittedly a silly choice, will change )

    movingPoints = [pointData(1).xyPoints;pointData(1).middleGrid];
    fixedPoints = [pointData(2).xyPoints;pointData(2).middleGrid];
    
    %=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=
    
    midGrid_m{cornerN} = movingPoints(end-2:end,:);
    midGrid_f{cornerN} = fixedPoints(end-2:end,:);
       
    evenlydistMoving(counter:counter+pointsPerSide-1,:) = movingPoints(1:pointsPerSide,:);
    evenlydistFixed(counter:counter+pointsPerSide-1,:) = fixedPoints(1:pointsPerSide,:);
    
    counter = counter+pointsPerSide;
end
%% ||~~~~-~~~~||~~ END LOOP ~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~

%% Now use middle Grid points to locate 9 internal control points

mside1 = [midGrid_m{1};midGrid_m{2}];
mside2 = [midGrid_m{3};midGrid_m{4}];
gridpoints_m = [mside1, mside2];

fside1 = [midGrid_f{1}; midGrid_f{2}];
fside2 = [midGrid_f{3}; midGrid_f{4}];
gridpoints_f = [fside1,fside2];

moving_coordinates = solveForGridPoints(gridpoints_m);
fixed_coordinates = solveForGridPoints(gridpoints_f);

cp_moving = [evenlydistMoving; moving_coordinates];
cp_fixed = [evenlydistFixed; fixed_coordinates];

%% Refine control points placement as necessary using this GUI!
close all force
[cp_moving, cp_fixed] = visualizeControlPoints_andResetManually(cp_moving, cp_fixed, MOVING_gray, IMG_gray);

close all 
figure;
showMatchedFeatures(MOVING_gray, IMG_gray, cp_moving, cp_fixed);

%% NOW WE USE CPOINTS TO DEFINE 3x Geometric nonrigid transformations

%this is allegedly a great preprocceesing step for multimodal registration
MOVING_gray = imhistmatch(MOVING_gray, IMG_gray);

nP = round(length(cp_moving)*0.9); %number of points to include in the local weighted means
if nP>=7
    tform1_lwn = cp2tform(cp_moving, cp_fixed, 'lwm', nP);
else
    nP = length(cp_moving)-1; %i.e. just shy of 100% my reasoning being even if you have zero edge points, 9 interior points -1 is still > 8 7
    tform1_lwn = cp2tform(cp_moving, cp_fixed, 'lwm', nP);
end
    tform2_poly2 = cp2tform(cp_moving, cp_fixed, 'polynomial', 2);
    tform3_poly3 = cp2tform(cp_moving, cp_fixed, 'polynomial', 3);

imReg1 = imtransform(MOVING_gray,tform1_lwn,'Xdata',[1 size(MOVING_mask,2)],'YData',[1 size(MOVING_mask,1)],'XYscale',1, 'FillValue', 1);
imReg2 = imtransform(MOVING_gray,tform2_poly2,'Xdata',[1 size(MOVING_mask,2)],'YData',[1 size(MOVING_mask,1)],'XYscale',1,'FillValue', 1);
imReg3 = imtransform(MOVING_gray,tform3_poly3,'Xdata',[1 size(MOVING_mask,2)],'YData',[1 size(MOVING_mask,1)],'XYscale',1,'FillValue', 1);

% call GUI to select best nonrigid transformation
choice = evaluate3nonRigidTransformations(MOVING_gray, imReg1, imReg2, imReg3, IMG_gray);
disp(num2str(choice));
assert(isscalar(choice));
switch choice
    case 11
        nearlyRegisteredMovingImage = imReg1;
        tform = tform1_lwn;
    case 22
        nearlyRegisteredMovingImage = imReg2;
        tform = tform2_poly2;
    case 33
        nearlyRegisteredMovingImage = imReg3;
        tform = tform3_poly3;
end

MOVINGMaskReg = imtransform(MOVING_mask, tform, 'Xdata',[1 size(MOVING_mask,2)],'YData',[1 size(MOVING_mask,1)],'XYscale',1, 'FillValue', 0);

%turn the warning back on you turned off at the beginning of the function
warning('on', id);

%close anything still open and be sure theyre really closed
close all force
pause(0.5);

%% FINAL REGISTRATION STEP!! Diffeomorphic demons
tic;
[D, movedIm_partway] = imregdemons(nearlyRegisteredMovingImage, IMG_gray, [750, 500, 30], 'PyramidLevels', 3, 'DisplayWaitbar', false, 'AccumulatedFieldSmoothing' , 1.25);
[D2, ~] = imregdemons(IMG_gray, movedIm_partway, [750, 500, 30], 'PyramidLevels', 3, 'DisplayWaitbar', false, 'AccumulatedFieldSmoothing' , 1.25);
toc;
D2inv = D2.*-1;
movedIm = imwarp(movedIm_partway, D2inv);
movedMask_partway = imwarp(MOVINGMaskReg, D);
movedMask = imwarp(movedMask_partway, D2inv);

% sz_im = size(movedIm, 1:2);
% [x,y ] = meshgrid(1:1:sz_im(2), 1:1:sz_im(1));
% figure;
% quiver(x,y,D2(:,:,1), D2(:,:,2));

%sum of squared differences
ssd0 = immse(nearlyRegisteredMovingImage, IMG_gray) * numel(nearlyRegisteredMovingImage);
deltaSSD = ssd0 - immse(movedIm, IMG_gray) * numel(movedIm);
partway_deltaSSD = ssd0 - immse(movedIm_partway, IMG_gray) * numel(movedIm_partway);

% Pearsons correlation coefficient
cor0 = corr2(nearlyRegisteredMovingImage, IMG_gray);
deltaCor = corr2(movedIm, IMG_gray)/cor0;
partway_deltaCor = corr2(movedIm_partway, IMG_gray)/cor0;

ccMoved= corr2(movedIm,nearlyRegisteredMovingImage);
ccFinalIms = corr2(movedIm, IMG_gray);
ccFinalMasks= corr2(movedMask, IMG_mask);

%Jaccard Distance: intersection over union
jacc_mask_partway = 1 - sum(sum(movedMask_partway & IMG_mask)) / sum(sum(movedMask_partway | IMG_mask));
jacc_masks = 1 - sum(sum(movedMask & IMG_mask)) / sum(sum(movedMask | IMG_mask));

disp(strcat(' correlation between final moving and fixed images ', num2str(ccFinalIms), 'and between their respective masks: ', num2str(ccFinalMasks)));
disp(strcat(' Pearsons correlation between moved images, before and after both demons is : ', num2str(ccMoved)));
disp(strcat(num2str(deltaCor), ' is the percent change in Pearsons Corr. Coef before and after both Demons (=after / before)'));
disp(strcat(num2str(deltaSSD), ' is the change in sum squared differences before and after both Demons (=before-after)'));
disp(strcat(num2str(partway_deltaCor), ' is the percent change in Pearsons Corr. Coef before and after just the 1st Demons (=after / before)'));
disp(strcat(num2str(partway_deltaSSD), ' is the change in sum squared differences before and after 1st Demons (=before-after)'));
disp(strcat(num2str(jacc_masks), ' is the jaccard distance of the image masks after Both runs'));
disp(strcat(num2str(jacc_masks), ' is the jaccard distance of the image masks after Both runs'));
disp(strcat(num2str(jacc_mask_partway), ' is the jaccard distance of the image masks PARTWAYS'));




%% Visualization #3 - does not work yet!!!

% create all of your visualizations -- lots of fused overlays
movingIm_fuse = imfuse(movedIm, nearlyRegisteredMovingImage, 'falsecolor', 'ColorChannels', [2,2,1],  'scaling', 'none');
partwayImg_fuse = imfuse(movedIm_partway, IMG_gray, 'falsecolor', 'ColorChannels', [2,1,2],'scaling', 'none');
finalImg_fuse = imfuse(movedIm, IMG_gray, 'falsecolor', 'ColorChannels', [2,1,2], 'scaling', 'none');

movingMasks_fuse = imfuse(movedMask, MOVINGMaskReg, 'falsecolor', 'ColorChannels', [2,2,1], 'scaling', 'none');
partwayMask_fuse = imfuse(movedMask_partway, IMG_mask, 'falsecolor', 'ColorChannels', [2,1,2], 'scaling', 'none');
finalMask_fuse = imfuse(movedMask, IMG_mask, 'falsecolor', 'ColorChannels', [2,1,2],'scaling', 'none');

MOVING_startPart2_Im = imfuse(movedIm, myNEWimages{2},  'falsecolor', 'ColorChannels', [2,1,2], 'scaling', 'none');
MOVING_startPart2_mask = imfuse(movedMask, myNEWimages{1}, 'falsecolor', 'ColorChannels', [2,1,2], 'scaling', 'none');

Dfield = imresize(imbinarize(checkerboard(100)), size(D, 1:2));
Dfield_part1 = padarray(imwarp(Dfield, D), [10, 10], 1, 'both');
Dfield_part2alone = padarray(imwarp(Dfield,D2inv), [10, 10], 1, 'both');
Dfield_combined = padarray(imwarp(Dfield_part1, D2inv), [10, 10], 1, 'both');

% now create the figure;

f3 = uifigure('Visible', 'on');
set(f3, 'units', 'pixels');

%make the GUI bigger
%f3_position = get(f3, 'Position');
f3_position = [10, 10, 1350 800];
set(f3, 'Position', f3_position);

gl_3 = uigridlayout(f3, [3,5], 'RowHeight', {'1x', '1.1x', 50});


butClose = uibutton(gl_3,'push', ...
    'Text','Close The Visualization?',...
    'ButtonPushedFcn', @(~,~) butCloseFcn);
butClose.Layout.Row = 3; butClose.Layout.Column = 1;

axMovingImsDemons = uiaxes(gl_3);  title(axMovingImsDemons, 'Moving Ims before/after demons');
axPartwaysDemons = uiaxes(gl_3); title(axPartwaysDemons, {'Fixed Im and Moving Im';'after first half of demons'});
axFinalDemons = uiaxes(gl_3); title(axFinalDemons, {'Fixed Im and Moving Im';'after both demons run'});
axAllPart2 = uiaxes(gl_3); title(axAllPart2, 'Comparison input to Part 2 Reg and Output')
axDfields= uiaxes(gl_3); title(axDfields,{ 'Visualize the Demons field D.' ;'BotRight is the combined'})
%axDfields.PositionConstraint = 'outerposition';

 axMovingImsDemons.Layout.Column = [3, 4];   axMovingImsDemons.Layout.Row = 1;
 axPartwaysDemons.Layout.Column = [1, 2] ;   axPartwaysDemons.Layout.Row = 2; 
 axFinalDemons.Layout.Column = [3, 4]    ;   axFinalDemons.Layout.Row = 2;
 axAllPart2.Layout.Column = [1, 2]       ;   axAllPart2.Layout.Row = 1;
 axDfields.Layout.Column = 5             ;   axDfields.Layout.Row = [1, 3];
 
 axis image;

 
imshowpair(MOVING_startPart2_mask, MOVING_startPart2_Im, 'montage', 'Parent', axAllPart2);
imshowpair(movingMasks_fuse, movingIm_fuse, 'montage', 'Parent', axMovingImsDemons);
imshowpair(partwayMask_fuse, partwayImg_fuse, 'montage', 'Parent', axPartwaysDemons);
imshowpair(finalMask_fuse, finalImg_fuse, 'montage', 'Parent', axFinalDemons);
montage({Dfield_part1, Dfield_part2alone, Dfield_combined}, 'Size', [3,1], 'Parent', axDfields);

drawnow nocallbacks;
f3.Visible = 'on';

uiwait
close all force

% M_Im and moved Mask are output variables for the function = they get sent
% back to registration part 1 where they are saved!!
end




function choice = evaluate3nonRigidTransformations(MOVING_gray, imReg1, imReg2, imReg3, IMG_gray)

hFig = uifigure(...
    'Name', 'Registration', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'Toolbar', 'none',...
    'Visible', 'off');

set(hFig, 'units', 'pixels');
pos = get(hFig, 'Position');
pos(3:4) = [800 600];
set(hFig, 'Position', pos);

%obj.Handles.Figure = hFig;
gl = uigridlayout(hFig,[3, 3],...
    'RowHeight', {'1x',40 , 40});

ax1 = uiaxes('Parent', gl);
ax1.Layout.Row = 1;
ax1.Layout.Column = 1;

ax2 = uiaxes('Parent', gl);
ax2.Layout.Row = 1;
ax2.Layout.Column = 2;

ax3 = uiaxes('Parent', gl);
ax3.Layout.Row = 1;
ax3.Layout.Column = 3;

im1 = imfuse(imReg1, IMG_gray, 'falsecolor', 'scaling', 'none', 'ColorChannels', [1,2,2]);
im2 = imfuse(imReg2, IMG_gray, 'falsecolor', 'scaling', 'none', 'ColorChannels', [1,2,2]);
im3 = imfuse(imReg3, IMG_gray, 'falsecolor', 'scaling', 'none', 'ColorChannels', [1,2,2]);

imshow(im1,'Parent', ax1, 'Border', 'tight');
imshow(im2,'Parent', ax2, 'Border', 'tight');
imshow(im3,'Parent', ax3, 'Border', 'tight');

title(ax1, 'LWM method');
title(ax2, 'Polynomial, deg 2');
title(ax3, 'Polynomial, deg 3');

% sum of squared differences
ssd0 = immse(MOVING_gray, IMG_gray) * numel(MOVING_gray);
ssd1 = ssd0 - immse(imReg1, IMG_gray) * numel(imReg1);
ssd2 = ssd0 - immse(imReg2, IMG_gray) * numel(imReg2);
ssd3 = ssd0 - immse(imReg3, IMG_gray) * numel(imReg3);

% correlation
cor0 = corr2(MOVING_gray, IMG_gray);
cor1 = corr2(imReg1, IMG_gray)/cor0;
cor2 = corr2(imReg2, IMG_gray)/cor0;
cor3 = corr2(imReg3, IMG_gray)/cor0;

disp('better registration should result in >> 1 delta corr, less than 1 means image is now less aligned (=corr2(im_reg, fixde)/corr2(im, fxed))');
disp('better reg. should also result in higher change sum of squared differences (=SSD(im, fxed)-SSD(im_reg, fxed))');

text1 = strcat('deltaSSD = ' ,num2str(ssd1),'delta corr = ' ,num2str(cor1));
text2 = strcat('deltaSSD = ' ,num2str(ssd2),'delta corr = ' ,num2str(cor2));
text3 = strcat('deltaSSD = ' ,num2str(ssd3),'delta corr = ' ,num2str(cor3));

lbl1 = uilabel(gl, 'Text', text1);
lbl1.Layout.Column = 1;
lbl2 = uilabel(gl, 'Text', text2);
lbl2.Layout.Column = 2;
lbl3 = uilabel(gl, 'Text', text3);
lbl3.Layout.Column = 3;

lbl1.Layout.Row = 2; lbl2.Layout.Row = 2; lbl3.Layout.Row = 2;

btn1 = uibutton(gl, 'push', 'Text', 'LWM', 'ButtonPushedFcn', @button1Callback);
btn2 = uibutton(gl, 'push', 'Text', 'polynomial, deg2','ButtonPushedFcn', @button2Callback);
btn3 = uibutton(gl, 'push', 'Text', 'polynomial, deg3','ButtonPushedFcn', @button3Callback);

btn1.Layout.Row = 3;
btn1.Layout.Column = 1;
btn2.Layout.Row = 3;
btn2.Layout.Column = 2;
btn3.Layout.Row = 3;
btn3.Layout.Column = 3;

set(hFig, 'Visible', 'on');

uiwait;
% retrieve app data holding user selection from GUI
choice = getappdata(0, 'mySelection');

%reset app data for next time!
setappdata(0, 'mySelection', []);
end

function button1Callback(~, ~)
setappdata(0, 'mySelection', 11);

uiresume;

end

function button2Callback(~, ~)
setappdata(0, 'mySelection', 22);

uiresume;

end

function button3Callback(~, ~)
setappdata(0, 'mySelection', 33);

uiresume;

end

function butCloseFcn(~,~)
uiresume;
end

function [cp_moving, cp_fixed] = visualizeControlPoints_andResetManually(cp_moving, cp_fixed, MOVING_gray, IMG_gray)

close all force

ff = figure;
imshowpair(MOVING_gray, IMG_gray, 'blend');
% Define the color order for the plots
colors = {'r', 'g'};
cp = [cp_moving; cp_fixed];
sz = size(cp_moving, 1); %this necessarily is the same value as size(cp_fixed, 1)
% Plot the control points on both images
hold on;
plot(cp(1:sz,1), cp(1:sz,2), [colors{1} '*'], 'MarkerSize', 10);
plot(cp(sz+1:2*sz,1), cp(sz+1:2*sz,2), [colors{2} '*'], 'MarkerSize', 10);
hold off;

disp('Please click on any misaligned control points to adjust their position. Press any key besides enter to finish.');
loopy=1;
while loopy == 1
    [x, y, butt] = ginput(1);
    if isempty(butt) || any(~isempty(butt)& butt~=1) || numel([x,y])==0  % exit loop if button other than left-click is pressed
        loopy = 9;
        continue
    else
        
        % Find the nearest control point to the clicked position
        try
            distances = sqrt(sum(bsxfun(@minus, [x y], cp).^2, 2));
        catch
             loopy = 9;
            continue
        end
        [~, idx] = min(distances);
        
        
        % Determine which set of control points the selected point belongs to
        if idx <= sz
            curr_cp = cp(1: sz, 1:2);
            flag = 1;
            %  curr_h = h_moving;
        elseif idx > sz
            curr_cp = cp(sz+1:end, 1:2);
            flag = 2;
            %curr_h = h_fixed;
            idx = idx - sz;
        end
        oldx = curr_cp(idx, 1);
        oldy = curr_cp(idx, 2);
               
        cla;
        imshowpair(MOVING_gray, IMG_gray, 'blend');
        %instead of deleting the point and messing up our index values. we just
        %move it out of frame,
        if flag == 1
            cp(idx, 1:2) = [-10, -10];
        else
            cp(idx+sz, 1:2)= [-10, -10];
        end
        
        %replot without that 1 point
        hold on
        plot(cp(1:sz,1), cp(1:sz,2), [colors{1} '*'], 'MarkerSize', 10);
        plot(cp(sz+1:2*sz,1), cp(sz+1:2*sz,2), [colors{2} '*'], 'MarkerSize', 10);
        %replace the old point a blue marker as filler
        hOld = plot(oldx, oldy, 'bo', 'MarkerSize', 10);
        hold off
        
        % Prompt the user to adjust the control point's position
        new_pos = ginput(1);
        ff.Visible = 'off';
        [ii, ~] = find(cp<0, 1, 'first'); % find the first negative value of cp, whould should be your point
        %assin the new value of the replaced coordinate into cp
        cp(ii, :) = new_pos(1:2);
        
        delete(hOld);
        cla; imshowpair(MOVING_gray, IMG_gray, 'blend');
        hold on;
        plot(cp(1:sz,1), cp(1:sz,2), [colors{1} '*'], 'MarkerSize', 10);
        plot(cp(sz+1:2*sz,1), cp(sz+1:2*sz,2), [colors{2} '*'], 'MarkerSize', 10);
        hold off;
        ff.Visible = 'on';
    end
    clear x y
end
cp_moving = cp(1:size(cp_moving, 1),1:2);
cp_fixed = cp(size(cp_moving, 1)+1:end,1:2);
end



