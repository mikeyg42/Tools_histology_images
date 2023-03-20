function [movingPoints, fixedPoints, MOVINGmask_fin, IMGmask_fin] = selectCornersGUI(MOVINGmask, IMGmask, MOVINGgray,IMGgray)
   
%% NOTE: mask #1 is MOVING
%% NOTE: mask #2 is FIXED_IMG

% Step 1: get rid of all but the ONE largest blob in the mask image
unadulteratedMasks = {MOVINGmask, IMGmask};
MyMasks= cell(1:2);
for jj = 1:2
    CurrentMask = unadulteratedMasks{jj};
    CC_MOVING = bwconncomp(CurrentMask);
    nPixs_MOVING = cellfun(@numel,CC_MOVING.PixelIdxList);
    [~,idx] = max(nPixs_MOVING);
    rng = 1:1:size(nPixs_MOVING, 2);
    differentIndex = setdiff(rng, idx);
    
    for pp = 1:numel(differentIndex)
        idxx = CC_MOVING.PixelIdxList{differentIndex(pp)};
        CurrentMask(idxx) = 0;
    end    
    MyMasks{jj}=imfill(CurrentMask, 4, 'holes');
end

% Step 2: convert mask into a polygon, and then use the nearestVertex
% function to query the 4 corners of the whole mask to locate which vertecies
% of the polygon are closest to the corners. These vertecies are likely
% corners of tissue. 

[simpHulls, cornerPointsInfo, linesToCornersInfo] = putativeCorners(MyMasks);

close all force

fPoints = struct([]); %initialize
ys = cell(1,2);
xs = cell(1,2);
for jj = 1:2
    %% This calls the actual GUI! 
    results = setMaskCornerPoints(linesToCornersInfo(jj), MyMasks{jj});
    newFields = {results.topLeft, results.topRight, results.botRight, results.botLeft};
    
%In case the selected corner point is just outside of the mask, or has subpixel accuracy
    points = simpHulls{jj};
    for n = 1:4
        coords = newFields{n};
        newCorner = coords(2,:);
        [kk,~] = dsearchn(points,newCorner);
        fPoints(jj).coordinates{n} = points(kk,:);
    end
    
    mIdx = cell2mat(fPoints(jj).coordinates);
    ys{jj} = [mIdx(2); mIdx(4); mIdx(6); mIdx(8)];
    xs{jj} = [mIdx(1); mIdx(3); mIdx(5); mIdx(7)];
end

% marginally better w/ these 2 split for-loops bc 2nd GUI can immediately load after 1st is closed!!
for jj = 1:2
% adding in centroid information

    ss = regionprops(MyMasks{jj}, 'Centroid');
    centroid =[ss.Centroid];
    
    C_x = (centroid(1));
    C_y = (centroid(2));
    cornerPointsInfo(jj).centroid = [C_x, C_y];
 
end

fixedPoints = [xs{2}, ys{2}];
MP = [xs{1}, ys{1}];
% use cross correlation to fix the selected points (moves them max 4pi)
movingPoints = cpcorr(MP, fixedPoints, MOVINGgray, IMGgray);

MOVINGmask_fin = MyMasks{1};
IMGmask_fin = MyMasks{2};

end


function [simpHulls, cornerPointsInfo, linesToCornersInfo] = putativeCorners(masks)
%masks is a cell array of your masks. This function is called out
%SelectCornersGUI_wrapper.m

linesToCornersInfo = struct('topLeft', [],'topRight',[],'botRight', [],'botLeft', [], 'cornerDistances', []);
cornerPointsInfo= struct('topLeft', [], 'topRight', [], 'botRight', [], 'botLeft', []);
simpHulls = cell(1,2);
for ii = 1:numel(masks)
    maskIMG = masks{ii};
    mask1 = bwareafilt(maskIMG,1);
    [y, x] = find(mask1==1); % results in a list of all white pixels
    
    % find the outermost of your points
    idx = convhull(x, y);
    bx = x(idx);
    by = y(idx);
    % make points into a polyshape
    hull_pgon = polyshape(bx,by,'Simplify',true, 'KeepCollinearPoints' , true);
    simpHulls{ii} = hull_pgon.Vertices;        % you will need this cell array later on!
    [xb, yb] = boundary(hull_pgon);
    
    % simplify Hull even more using the REDUCEPOLY function
    simplep = reducepoly([xb, yb]);
    simplepgon = polyshape(simplep);
    coordsSimpleHull = simplepgon.Vertices;
    
    % we are looping through each point of our shape to see which vert is
    % closest to the corners
    sz_curIm = size(maskIMG, 1:2);
    cornerDist = struct('tLeft', [], 'tRight', [], 'bRight', [],'bLeft', []);
    for t = 1:size(coordsSimpleHull, 1)
        mIdx = abs(coordsSimpleHull(t, 1)-1);
        fIdx = abs(coordsSimpleHull(t, 2)-1);
        c = abs(coordsSimpleHull(t, 1)-sz_curIm(2));
        d = abs(coordsSimpleHull(t, 2)-sz_curIm(1));
        cornerDist(t).tLeft = hypot(mIdx,fIdx);
        cornerDist(t).bLeft = hypot(mIdx,d);
        cornerDist(t).tRight = hypot(c,fIdx);
        cornerDist(t).bRight  = hypot(c,d);
    end
    linesToCornersInfo(ii).cornerDistances = cornerDist;
    
    [~,tL] = min([cornerDist.tLeft]);
    [~,bL] = min([cornerDist.bLeft]);
    [~,bR] = min([cornerDist.bRight]);
    [~,tR] = min([cornerDist.tRight]);
    
    topLeftExtreme = [1, 1];
    topRightExtreme = [sz_curIm(2), 1];
    botRightExtreme = [sz_curIm(2), sz_curIm(1)];
    botLeftExtreme = [1, sz_curIm(1)];
    
    linesToCornersInfo(ii).topLeft = [topLeftExtreme; coordsSimpleHull(tL, :)];
    linesToCornersInfo(ii).topRight =  [topRightExtreme; coordsSimpleHull(tR, :)];
    linesToCornersInfo(ii).botRight = [botRightExtreme; coordsSimpleHull(bR, :)];
    linesToCornersInfo(ii).botLeft = [botLeftExtreme; coordsSimpleHull(bL, :)];
    
    cornerPointsInfo(ii).topLeft = coordsSimpleHull(tL, :);
    cornerPointsInfo(ii).topRight = coordsSimpleHull(tR, :);
    cornerPointsInfo(ii).botRight = coordsSimpleHull(bR, :);
    cornerPointsInfo(ii).botLeft = coordsSimpleHull(bL, :);
    
end
end



