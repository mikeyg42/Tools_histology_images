function [D, tform, M_Im, movedMask] = registerSerialSections_part2_nonrigid(varargin)

%% Part 2 of registration: nonrigid geometric transformation with control points placed EXCLUSIVELY programically .
% name of the game is placing control points programmically, and with great
% care.

%Using points, I try some nonrigid registration ( LWH / polynomials). After
%picking the best of those models, I do the final cajoling of the points
%together using a modified diffeomorphic demons approach

id = 'MATLAB:polyshape:repairedBySimplify';
warning('off', id);



%% load up variables saved 
if nargin ~= 3
    ld = load('/Users/jglendin/Dropbox - Michael/Dropbox/human ms code/my code/registration_data_3164_MCAM_PLP.mat',...
        '-mat');
    myNEWimages = ld.images_part1;
    fixedMat = ld.fixedMat;
    movingMat = ld.movingMat;
else
    myNEWimages = varargin{3};
    fixedMat = varargin{1};
    movingMat = varargin{2};
end

LESION_gray = myNEWimages{2};
LESION_mask = myNEWimages{1};
IMG_gray = myNEWimages{4};
IMG_mask = myNEWimages{3};

LESION_gray(~LESION_mask) = 1;
IMG_gray(~IMG_mask) = 1;

%% define evenly spaced points along each of the four edges -- will become the majority of control points

%----------------------------------------------------------------------------
% start first with the fixed image
%----------------------------------------------------------------------------
% 1. make hull
[ay, ax]=find(IMG_mask); 
hull_1 = convhull(ax, ay);

% 2. convert hull into polyshape
HULL_image1 = polybuffer(polyshape(ax(hull_1), ay(hull_1)),  0.1, 'JointType','miter');

% 3. get coordinaes of contour
img_Cont = contourc(im2double(IMG_mask), [0.5, 0.5]); 
img_Cont(:, 1) = [];
cxp = img_Cont(1,2:end)'; 
cyp = img_Cont(2,2:end)';
cxp(cxp>size(IMG_mask, 2)) = size(IMG_mask, 2);
cyp(cyp>size(IMG_mask, 1)) = size(IMG_mask, 1);
xdiff_m = [diff(cxp); cxp(end)-cxp(1)]; ydiff_m = [diff(cyp); cyp(end)-cyp(1)];
midx = ( xdiff_m>1 | ydiff_m>1);

% 4. query polyshape with the contour lines' coordinates
TFin1 = isinterior(HULL_image1, cxp, cyp);
cxp(~TFin1 | midx) = []; 
cyp(~TFin1 | midx) = []; 

% 5. use boundary and bwtraceboundary together to ensure points are in right order 
k2 = boundary(cxp, cyp, 0.9);
imgBW=poly2mask(cxp(k2), cyp(k2), size(IMG_mask, 1), size(IMG_mask,2));
bw3 = bwperim(imgBW);
[r0w,c0l] = find(bw3);
[~, jj] = pdist2([c0l,r0w],[1, 1], 'euclidean', 'Smallest', 1);
C = bwtraceboundary(imgBW, [r0w(jj),c0l(jj)], 'NE'); 

% 6. make the shapefrom the contour vertices that were not outliers!
pgon_fixed = polyshape(C(:,2), C(:,1),  'SolidBoundaryOrientation', 'cw');
pgon_fixed = simplify(pgon_fixed,'KeepCollinearPoints',1);
pgonFixed_simpler = rmholes(pgon_fixed); 

% +++++++++
% 1. make hull
[by, bx]=find(LESION_mask); 
hull_2 = convhull(bx, by);

% 2. convert hull into polyshape
HULL_image2 = polybuffer(polyshape(bx(hull_2), by(hull_2)), 0.1, 'JointType','miter'); 

% 3. get coordinaes of contour
lesion_Cont = contourc(im2double(LESION_mask), [0.5, 0.5]); 
lesion_Cont(:, 1) = [];
dxp = lesion_Cont(1,2:end)'; 
dyp = lesion_Cont(2,2:end)';
dxp(dxp>size(LESION_mask, 2)) = size(LESION_mask, 2);
dyp(dyp>size(LESION_mask, 1)) = size(LESION_mask, 1);
xdiff = [diff(dxp); dxp(end)-dxp(1)]; ydiff= [diff(dyp); dyp(end)-dyp(1)];
idx = ( xdiff>1 | ydiff>1);

% 4. query polyshape with the contour lines' coordinates
TFin2 = isinterior(HULL_image2, dxp, dyp);
dxp(~TFin2 | idx) = []; 
dyp(~TFin2 | idx) = []; 

% 5. use boundary and bwtraceboundary together to ensure points are in right order 
k = boundary(dxp, dyp, 0.9);
lesBW=poly2mask(dxp(k), dyp(k), size(LESION_mask, 1), size(LESION_mask,2));
bw2 = bwperim(lesBW);
[rw,cl] = find(bw2);
[~, ii] = pdist2([cl,rw],[1, 1],'euclidean', 'Smallest', 1);
B = bwtraceboundary(lesBW, [rw(ii),cl(ii)], 'NE'); 

% 6. make the shapefrom the contour vertices that were not outliers!
pgon_move = polyshape(B(:,2), B(:,1), 'SolidBoundaryOrientation', 'cw');
pgon_move = simplify(pgon_move,'KeepCollinearPoints',1);
pgonMOVING_simpler = rmholes(pgon_move);

%----------------------------------------------------------------------------
%% these polygons are not quite identical to mask. Therefore we need to adjust the corner points
% we want our corners back because we will use them considerably more

cornerIndxMoving = nearestvertex(pgonMOVING_simpler, movingMat(1:4,:));
cornerIndxFixed = nearestvertex(pgonFixed_simpler, fixedMat(1:4,:));

movingM_pgon = pgonMOVING_simpler.Vertices(cornerIndxMoving(1:4, 1), :);
fixedM_pgon = pgonFixed_simpler.Vertices(cornerIndxFixed(1:4, 1), :);

% if the diffference between matrix coordinate and the nearest point on
% plygon is more than 4 pixels, then CPCORR function won't fix it.
diffCorn = movingM_pgon - movingMat(1:4,:);
diffCorn_2 = fixedM_pgon - fixedMat(1:4,:); 
Midx = (diffCorn(:,1).^2+diffCorn(:,2).^2).^0.5 > 4;
Fidx = (diffCorn_2(:,1).^2+diffCorn_2(:,2).^2).^0.5 > 4;
for L = 1:4
    if Midx(L)==1
        movingM_pgon(L,:) = (movingM_pgon(L,:) + movingMat(L,:))./2;
    end
    if Fidx(L)==1
        fixedM_pgon(L,:) = (fixedM_pgon(L,:) + fixedMat(L,:))./2;
    end
end

fixedM = fixedM_pgon;
movingM = movingM_pgon;

% wiggle points so that they register better
movingM = cpcorr(movingM,fixedM,LESION_gray,IMG_gray); 
fixedM = cpcorr(fixedM,movingM, IMG_gray,LESION_gray);

% ----------- quick sanity check -----------
sumMoving = sum(movingM,2);
sumFixed = sum(fixedM, 2);
[~, id1] = min(sumMoving); [~, id3] = max(sumMoving);
if id1~=1 || id3 ~=3 || movingM(4, 1)>movingM(2,1) || movingM(4, 2)<movingM(2,2)
disp('error of order: MOVING'); pause(1);
end

[~, fid1] = min(sumFixed); [~, fid3] = max(sumFixed);
if fid1~=1 || fid3 ~=3 || fixedM(4, 2)<fixedM(2,2) || fixedM(4, 1)>fixedM(2,1)
disp('error of order: FIXED'); pause(1);
end
% -----------  end sanity check  -----------

%% Use curve fitting to place evenly spaced points along each edge of the tissue


%% SET NUM POINTS ALONG AN EDGE
pointsPerSide = 16; % multiplied by 4 +9 gives total CP's
%%
%preallocate
evenlydistMoving = zeros(pointsPerSide*4,2, 'double'); 
evenlydistFixed = zeros(pointsPerSide*4,2,'double');
MovingPolynomialVals =zeros(pointsPerSide*4,2,'double'); 
FixedPolynomialVals = zeros(pointsPerSide*4,2,'double');
distancefromcurvetopolygon_m = zeros(pointsPerSide, 1,'double');
distancefromcurvetopolygon_f = zeros(pointsPerSide, 1,'double');

tic;

counter = 1;
for cornerN = 1:4
%corners are numbered 1 -> 4 CW starting in top left. edges are named by the corners they span. 1-2, 2-3, 3-4, and 4-1.    
    if cornerN~=4
        cornerN_plus1 = cornerN+1;
    else
        cornerN_plus1=1;
    end
% make sure corners are in the right order in moving image
    cornsMoving = [movingM(cornerN, 1:2); movingM(cornerN_plus1, 1:2)];
    starting = cornsMoving(2, :);
    ending = cornsMoving(1,:);
    
    st = nearestvertex(pgonMOVING_simpler, starting(1,1), starting(1, 2));
    ed = nearestvertex(pgonMOVING_simpler, ending(1,1), ending(1, 2));

    if ed > st
        PointsToFit = pgonMOVING_simpler.Vertices(st:1:ed, 1:2); 
    else %this should only happen once, 
        PointsToFit = [pgonMOVING_simpler.Vertices(st:end, 1:2); pgonMOVING_simpler.Vertices(1:ed, 1:2)];
    end
    

    % repeat with fixed
    cornsFixed = [fixedM(cornerN, 1:2); fixedM(cornerN_plus1,1:2)];
    startingF = cornsFixed(2, :);
    endingF = cornsFixed(1,:);

    st2 = nearestvertex(pgonFixed_simpler, startingF(1,1), startingF(1, 2));
    ed2 = nearestvertex(pgonFixed_simpler, endingF(1,1), endingF(1, 2));
    if ed2 > st2
        PointsToFit_fix = pgonFixed_simpler.Vertices(st2:ed2, 1:2); 
    else %this should only happen once, when the highest index is reached and the numbers restart at 1
        PointsToFit_fix = [pgonFixed_simpler.Vertices(st2:end, 1:2) ;pgonFixed_simpler.Vertices(1:ed2, 1:2)];
    end

    
    %% call the curve fitting script

    pointData  = curveFittingOfTissueBorders(pointsPerSide, cornerN, movingM, fixedM, PointsToFit, PointsToFit_fix, IMG_gray, LESION_gray);
    
%%
    xy_moving= pointData(1).xyPoints; 
    xy_fixed = pointData(2).xyPoints;
    
    MovingPolynomialVals(counter:counter+pointsPerSide-1,:) = xy_moving(1:end-1,:);
    FixedPolynomialVals(counter:counter+pointsPerSide-1,:) = xy_fixed(1:end-1,:);   
    
%reformat middleGridPoints, which will be at most 9 potential internal control points
% two consecutive sides will need to be reversed order. it is arbitrary which 2
    if cornerN == 3 || cornerN ==4
        midGrid_m{cornerN} = flipud(pointData(1).middleGrid);
        midGrid_f{cornerN} = flipud(pointData(2).middleGrid);
    else
        midGrid_m{cornerN} = pointData(1).middleGrid;
        midGrid_f{cornerN} = pointData(2).middleGrid;
    end

%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=
movingContour = [B(:,2), B(:,1)];
fixedContour =  [C(:,2), C(:,1)];

    for jj = 1:pointsPerSide %first and last corner we just MAKE be the corners
        if jj == 1
            xy_moving(jj, 1:2) = starting;
            xy_fixed(jj, 1:2) = startingF;
            continue;
            
        elseif jj == pointsPerSide
            xy_moving(jj, 1:2) = ending; % the "last" point SHOULD ALWAYS be the next corner point! but we wait to get rid of it til end of loop!
            xy_fixed(jj, 1:2) = endingF; 
            continue;
            
        elseif mod(cornerN, 2) == 1  
            flagg = 1; 
            while flagg == 1
            indx_m = find(movingContour(:,1) == round(xy_moving(jj, 1)));
            
            nPts_m = numel(indx_m); % number of intersections of your xy_moving point with contour
                switch nPts_m
                    case 0
                        xy_moving(jj, 1) = round(xy_moving(jj, 1)+0.6); %shift the x coordinate half a pixel to the right, then try again
                        flagg = 1; %NOT DONE
                    case 1
                        distancefromcurvetopolygon_m(jj, 1) = pdist2(xy_moving(jj, 1:2), movingContour(indx_m, 1:2), 'euclidean');
                        xy_moving(jj, 2) = movingContour(indx_m, 2); %set the 
                        flagg=2; %DONE
                    otherwise % ie. " > 1 " THIS IS MOST CASES!!!!
                        hits = movingContour(indx_m,1:2);
                        [distances, Id_multi] = pdist2(hits, xy_moving(jj,1:2), 'euclidean', 'Smallest', 1);
                        distancefromcurvetopolygon_m(jj, 1) = distances(1, 1);
                        xy_moving(jj, 2) = hits(Id_multi,2);
                        flagg = 2; %DONEred
                end
            end
% ^ moving            
%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=
% v fixed

            flage = 1;
            while flage ==1
                indx_f = find(fixedContour(:,1) == round(xy_fixed(jj, 1)));
                nPts_f = numel(indx_f);
                switch nPts_f
                    case 0
                        xy_fixed(jj, 1) = round(xy_fixed(jj, 1)+0.6); 
                        flage = 1; %NOT DONE
                    case 1
                        distancefromcurvetopolygon_f(jj, 1) = pdist2(xy_fixed(jj, 1:2), fixedContour(indx_f, 1:2), 'euclidean');
                        xy_fixed(jj, 2) = fixedContour(indx_f,2);
                        flage = 2; %DONE
                    otherwise % ie. " > 1 "
                        hitts = fixedContour(indx_f,1:2);
                        [distys, Id2_multi] = pdist2(hitts, xy_fixed(jj, 1:2), 'euclidean', 'Smallest', 1);
                        distancefromcurvetopolygon_f(jj, 1) = distys(1, 1);
                        xy_fixed(jj, 2) = hitts(Id2_multi,2);
                        flage = 2; %DONE
                end
            end
  %                                            ^ ^
  % -=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-|%|-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=
  % +-=%=-+-=%=-+-=% this is identical to the above,with x and y coordinates swapped %=-+-=%=-+-
  % -+-=%=-+-=%=-+-=%|-|-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+
  %                  V V          
        elseif mod(cornerN, 2) == 0 %ie QQ = 2 or 4  
            
            flaggy =1;
            while flaggy ==1
            indx_f = find(fixedContour(:,2) == round(xy_fixed(jj, 2)));
            nPts_f = numel(indx_f);
                switch nPts_f
                    case 0
                        xy_fixed(jj, 2) = round(xy_fixed(jj, 2)+0.6); 
                        flaggy = 1; %NOT DONE
                    case 1
                        distancefromcurvetopolygon_f(jj, 1) = pdist2(xy_fixed(jj, 1:2), fixedContour(indx_f, 1:2), 'euclidean');
                        xy_fixed(jj, 1) = fixedContour(indx_f,1);
                        flaggy=2;
                    otherwise % ie. " > 1 "
                        hits = fixedContour(indx_f,1:2);
                        [dists, Id_multi] = pdist2(hits, xy_fixed(jj, 1:2), 'euclidean', 'Smallest', 1);
                        distancefromcurvetopolygon_f(jj, 1) = dists(1, 1);
                        xy_fixed(jj, 1) = hits(Id_multi,1);
                        flaggy=2;
                end
            end
% ^ fixed             
%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=
 % v moving
            flags = 1;
            while flags == 1
            indx_m = find(movingContour(:,2) == round(xy_moving(jj, 2)));
            nPts_m = numel(indx_m);
                switch nPts_m
                    case 0
                        xy_moving(jj, 2) = round(xy_moving(jj, 2)+0.6); 
                        flags = 1; %NOT DONE
                    case 1
                        distancefromcurvetopolygon_m(jj, 1) = pdist2(xy_moving(jj, 1:2), movingContour(indx_m, 1:2), 'euclidean');
                        xy_moving(jj, 1) = movingContour(indx_m,1);
                        flags=2;
                    otherwise % ie. " > 1 "
                        hits = movingContour(indx_m,1:2);
                        [distance, Id_multi] = pdist2(hits, xy_moving(jj, 1:2), 'euclidean', 'Smallest', 1);
                        distancefromcurvetopolygon_m(jj, 1) = distance(1, 1);
                        xy_moving(jj, 1) = hits(Id_multi,1);
                        flags=2;
                end
            end 
        end
 %=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=%=-+-=
    end %end looping through points
    
    %% visual #1+2 - shows once per side
    
%     figure(1); imshow(LESION_gray); hold on; plot(pgonMOVING_simpler); visboundaries(LESION_mask,'color', 'green', 'LineStyle', ':', 'LineWidth', 1, 'EnhanceVisibility', false); plot(PointsToFit(:,1), PointsToFit(:,2), 'r-', 'LineWidth', 3);
%    
%     figure; imshow(IMG_gray); hold on; plot(pgonFixed_simpler); visboundaries(IMG_mask,'color', 'green', 'LineStyle', ':', 'LineWidth', 1, 'EnhanceVisibility', false); plot(PointsToFit_fix(:,1), PointsToFit_fix(:,2), 'r-', 'LineWidth', 3);
%     drawnow;
%     
%     
%     
%     f2 = figure('Visible', 'off'); 
%     axy1 = {ax1, ax3, ax5, ax7};
%     axy2 = {ax2, ax4, ax6, ax8};
%     gl = griddedlayout(f2, 'Size', [4, 4]);
%     axy1{cornerN} = axes(gl, 'Layout.Row', cornerN, 'Layout.Column', 1);
%     
%     axy2{cornerN} = axes(gl, 'Layout.Row', cornerN, 'Layout.Column', 2);
% 
%     imshow(IMG_gray,  'Parent', axy2{cornerN}); hold on; plot(xy_fixed(1,1), xy_fixed(1, 2), 'm*',xy_fixed(2:end,1), xy_fixed(2:end,2), 'co');
%     
%      imshow(LESION_gray, 'Parent', axy3); 
%      hold on; 
%      plot(pgonMOVING_simpler);
%    
%      plot(xy_moving(1,1), xy_moving(1, 2), 'm*',xy_moving(:,1), xy_moving(:,2), 'co');
%    hold off
%    
%    
%     
%     f2.Visible = 'on';
%     drawnow limitrate nocallbacks
%     pause(5);
    
    %%
    fixed_residuals = sum(distancefromcurvetopolygon_f(:))/pointsPerSide;
    moving_residuals = sum(distancefromcurvetopolygon_m(:))/pointsPerSide;
   
    disp(strcat('FIXED: Avg. distance from point to curve: ', num2str(fixed_residuals)));
    disp(strcat('MOVING: Avg. distance from point to curve: ', num2str(moving_residuals)));

    evenlydistMoving(counter:counter+pointsPerSide-1,:) = xy_moving(1:end-1,:);
    evenlydistFixed(counter:counter+pointsPerSide-1,:) = xy_fixed(1:end-1,:);
    
    counter = counter+pointsPerSide;
end
%% ||~~~~-~~~~||~~ END HUGE LOOP ~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~~~~||~~~~-~

toc;

cPointsMoving = cpcorr(evenlydistMoving, evenlydistFixed, LESION_gray, IMG_gray);
cPointsFix = cpcorr(evenlydistFixed, cPointsMoving, IMG_gray, LESION_gray);

%% Now use middle Grid points to define the internal control points

mside1 = [midGrid_m{1};fliplr(midGrid_m{2})];
mside2 = [flipud(midGrid_m{3});flipud(fliplr(midGrid_m{4}))];
gridpoints_m = [mside1, mside2];

fside1 = [midGrid_f{1};fliplr(midGrid_f{2})];
fside2 = [flipud(midGrid_f{3});flipud(fliplr(midGrid_f{4}))];
gridpoints_f = [fside1,fside2];

moving_coordinates = solveForGridPoints(gridpoints_m);
fixed_coordinates = solveForGridPoints(gridpoints_f);

cp_moving = [cPointsMoving; moving_coordinates];
cp_fixed = [cPointsFix; fixed_coordinates];

%% DONE WITH CONTROL POINT PLACEMENT!!

figure; showMatchedFeatures(LESION_gray, IMG_gray, cp_moving, cp_fixed);

% IT WOULD BE COOL IF YOU COULD DELETE AN ERRANT POINT OR TWO PAIR FROM
% THIS IAMGE SOMEHOW?


%% NOW WE USE CPOINTS TO DEFINE 3x Geometric nonrigid transformations

%this is allegedly a great preprocceesing step for multimodal registration
LESION_gray = imhistmatch(LESION_gray, IMG_gray);

nP = round(length(cp_moving)*0.9); %number of points to include in the local weighted means
tform1_lwn = cp2tform(cp_moving, cp_fixed, 'lwm', nP);
tform2_poly2 = cp2tform(cp_moving, cp_fixed, 'polynomial', 2);
tform3_poly3 = cp2tform(cp_moving, cp_fixed, 'polynomial', 3);

imReg1 = imtransform(LESION_gray,tform1_lwn,'Xdata',[1 size(LESION_mask,2)],'YData',[1 size(LESION_mask,1)],'XYscale',1, 'FillValue', 1);
imReg2 = imtransform(LESION_gray,tform2_poly2,'Xdata',[1 size(LESION_mask,2)],'YData',[1 size(LESION_mask,1)],'XYscale',1,'FillValue', 1);
imReg3 = imtransform(LESION_gray,tform3_poly3,'Xdata',[1 size(LESION_mask,2)],'YData',[1 size(LESION_mask,1)],'XYscale',1,'FillValue', 1);

% call GUI to select best nonrigid transformation
choice = evaluate3nonRigidTransformations(imReg1, imReg2, imReg3, IMG_gray);

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

LesionMaskReg = imtransform(LESION_mask, tform, 'Xdata',[1 size(LESION_mask,2)],'YData',[1 size(LESION_mask,1)],'XYscale',1, 'FillValue', 0);

%turn the warning back on you turned off at the beginning of the function
warning('off', id);

%close anything still open and be sure theyre really closed 
close all force
pause(0.5);

%% FINAL REGISTRATION STEP!! Diffeomorphic demons 
[D, M_Im] = imregdemons(nearlyRegisteredMovingImage, IMG_gray, [48, 36, 12, 4], 'PyramidLevels', 4, 'DisplayWaitbar', false);
movedMask = imwarp(LesionMaskReg, D);

%% Visualization #3

f3 = uifigure('Visible', 'off');
gl_3 = uigridlayout(f3, [2,3]);
gl_3.RowHeight = {'1x',20};
butClose = uibutton(gl_3,'push', ...
    'Text','Close The Visualization?',...
    'ButtonPushedFcn', @(~,~) butCloseFcn);

butClose.Layout.Row = 2;butClose.Layout.Column = 2;

axImgs = uiaxes(gl_3);  title(axImgs, 'Fixed and Moving Images')
axMasks = uiaxes(gl_3); title(axMasks, 'Fixed and Moving MASKS')
axChange = uiaxes(gl_3); title(axChange, 'Moving masks before and after demons algo')

axChange.Layout.Column = 1; axImgs.Layout.Column = 3; axMasks.Layout.Column = 2;
axChange.Layout.Row = [1, 2]; axImgs.Layout.Row = [1, 2]; axMasks.Layout.Row = 1;

imshowpair(M_Im, IMG_gray, 'falsecolor', 'ColorChannels', [1,2,2], 'Parent', axImgs, 'scaling', 'none');
imshowpair(movedMask, fixedMask,'falsecolor', 'ColorChannels', [2,1,2], 'Parent', axMasks,  'none', 'scaling', 'none');
imshowpair(LesionMaskReg, movedMask, 'falsecolor', 'ColorChannels', [2,2,1], 'Parent', axChange, 'scaling', 'none');

f3.Visible = 'on';
uiwait
close all force

% M_Im and moved Mask are output variables for the function = they get sent
% back to registration part 1 where they are saved!!
end




function choice = evaluate3nonRigidTransformations(imReg1, imReg2, imReg3, IMG_gray)

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

lbl1 = uilabel(gl, 'Text', num2str(corr2(imReg1, IMG_gray)));
lbl1.Layout.Column = 1;
lbl2 = uilabel(gl, 'Text', num2str(corr2(imReg2, IMG_gray)));
lbl2.Layout.Column = 2;
lbl3 = uilabel(gl, 'Text', num2str(corr2(imReg3, IMG_gray)));
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
