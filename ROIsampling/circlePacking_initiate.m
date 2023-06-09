function [solution1, tracking, info] = circlePacking_initiate(binaryMask,ROIsize, offset, settings)
% Michael Glendinning, 2023

warning('off', 'all');

assert(isnumeric(offset)&isreal(offset))
if isempty(settings.chooseROI.roi_newRadiusBounds)
    % Parse the maximum radius/diameter
    maxDiameter = sqrt(2)*max(ROIsize);
    maxRadius = maxDiameter/2;
    minRadius = max(ROIsize)/2;
    minDiameter = max(ROIsize);
else
    maxRadius = settings.chooseROI.roi_newRadiusBounds(2);
    maxDiameter = maxRadius*2;
    minRadius = settings.chooseROI.roi_newRadiusBounds(1);
    minDiameter = minRadius*2;
end
info = struct( 'maxRadius', maxRadius, 'minRadius', minRadius, 'maxDiameter',maxDiameter, 'minDiameter', minDiameter);

% initiallize loop by precalculating initial guesses of the entirety of the hexagonal packing lattice starting locations
depth = ceil(size(binaryMask,1)/(maxRadius*sqrt(3)));
centers_xy = generateHexLattice(binaryMask, info.maxRadius-0.6, depth); % we need to subtract 0.6 so that we have >1 pixel of overlap
centers_xy = centers_xy + [offset, offset];
centers_xy = checkCircleCenterPoints(info.maxRadius, centers_xy, binaryMask);

counter = 1;

% This image we will update each time we add or change the position of a circle 
circleMaskAll_map = imcomplement(binaryMask);

% Calculate progress
totalBlobArea = numel(find(circleMaskAll_map(:)==0));
minCircArea = 1.571*max(ROIsize)^2/2;

maxNumCircles = floor(totalBlobArea/minCircArea);
tracking = zeros(maxNumCircles,4);

%% part 1 goal: arrange as many circles as will fit in a hexagonal lattice like manner
for j = 1:size(centers_xy, 1) 
    cir = images.roi.Circle('Radius', info.maxRadius,'Center', centers_xy(j, :));
    circleMask = createMask(cir, binaryMask);

    overlapTest = circleMask.* circleMaskAll_map;
    %I can allow a tiny bit of overlap just in case 
    if sum(overlapTest(:)) > 120
        continue
    end

    % use the ROI-generated binary image to update the circleMaskAll image
    circleMaskAll{counter} = circleMask;
    circleMaskAll_map = circleMaskAll_map | circleMask ;

    % Keep track of how much coverage you have so far
    tracking(counter,1) = numel(find(circleMaskAll_map(:)==0));

    % also store circle's center and radius
    tracking(counter, 2) = info.maxDiameter/2;
    tracking(counter, 3) = centers_xy(j, 1);
    tracking(counter, 4) = centers_xy(j, 2);

    counter = counter+1;
end
part1packed = counter-1;
disp(strcat('Part 1 packed : ' , num2str(part1packed), ' circles'));

% Initiallize the variables old and new coverage which define exit condition of loop
new_coverage = tracking(part1packed,1);
old_coverage = 0;
counter2 = 1;
while new_coverage~=old_coverage
    % replace old- with new and forget about old, now that you've started a new loop
    old_coverage = new_coverage;

    % use the image with white forground, circle overlays and background black to describe a
    % triangulation of remaining space in want of circles
    spaceRemaining = imcomplement(circleMaskAll_map);
    boundLeftOver = bwboundaries(spaceRemaining,  'CoordinateOrder' , 'xy');
    X = cellfun(@(x) x(:, 2), boundLeftOver, 'UniformOutput', false);
    Y = cellfun(@(x) x(:, 1), boundLeftOver, 'UniformOutput', false);
    ps = polyshape(Y,X, 'KeepCollinearPoints', true, 'Simplify', false);

    ps_regions = regions(ps);
    region_areas = area(ps_regions);
    ps_regions(region_areas<max(region_areas)) = []; % you only want the one largest boundary

    outside = rmholes(ps_regions);
    outpoints = outside.Vertices;
    outpoints = cleanedUpList(outpoints);

    inside = rmboundary(ps_regions,1);
    inpoints = inside.Vertices;
    inpoints = cleanedUpList(inpoints);

    verts_xy = [outpoints; inpoints];

    % constrains are described using edge matrix convention, i.e. first row = first edge,
    % defined by the two vertices that are part of that edge..s.t. it looks like [1, 2; 2, 3; 3 ... ect
    out = (1:size(outpoints,1))';
    in = (size(outpoints,1)+1:1:size(verts_xy,1))';
    constraints = vertcat([out, circshift(out,1,1)], [in, circshift(in,1,1)]);

    % Create delaunay triangulation, then use the isInterior object fcn to define triangulation
    % of only interior points! note how tr, although valuable, does not fullfill stringent def of a DT.
    dt = delaunayTriangulation(verts_xy, constraints);
    TF = isInterior(dt);
    tr = triangulation(dt.ConnectivityList(TF,:),dt.Points(:,1),dt.Points(:,2));

    % A circumcenter is the center of the only circle one can draw connecting all 3 vertices of a triangle.
    % (a fact that lies at the core of delauney triangulation)
    [circums, radii] = circumcenter(tr);
    Dcircum = pdist2(tracking(1:nnz(tracking)/4, 3:4), circums);

    goodradii = radii<= maxRadius & radii>=minRadius & min(Dcircum, [], 1)'> minRadius*2 & circums(: ,1)>1 & circums(: ,2)>1;
    if any(goodradii==1)
         radii(~goodradii) = [];
         circums(~goodradii, :) = [];

         [newradius, idx] = max(radii, [],'all', 'linear');
         centerPoint = circums(idx, :);
         cir = createMask(images.roi.Circle('Center',centerPoint, 'Radius', newradius), binaryMask);

         circleMaskAll_map = circleMaskAll_map | cir;
         
         new_coverage = numel(find(circleMaskAll_map(:)==0));
        if sum(cir(:))>pi*newradius^2*0.8
         circleMaskAll{part1packed+counter2} = cir;
         tracking(part1packed+counter2, 1) = new_coverage;
         tracking(part1packed+counter2, 2) = newradius;
         tracking(part1packed+counter2, 3) = centerPoint(:, 1);
         tracking(part1packed+counter2, 4) = centerPoint(:, 2);
         counter2 = counter2+1;
        end
    else
        break;
    end
end

% Convert raw pixel count to proportion 
tracking(:, 1) = tracking(:, 1)./totalBlobArea;

totalCircsPacked = numel(unique(round(tracking(:,1), 4)))-1;

part2packed = totalCircsPacked - part1packed;
disp(strcat('Part 2 packed : ' , num2str(part2packed), ' circles'));
finalCoveredArea = numel(find(circleMaskAll_map(:)==0));
disp(strcat('Total packed : ' , num2str(totalCircsPacked), ' with coverage of :', num2str(finalCoveredArea/totalBlobArea)));

solution1(1).all = circleMaskAll_map.*binaryMask;
solution1(1).individ = circleMaskAll;

warning('on', 'all');
end

function centers_xy = generateHexLattice(binaryMask, maxRadius, depth)
%recursively iterate row by row down the image, replicating and offsetting the rows such
%that they are tangent that they utilizing the most efficient stacking of 2D circles
%known. 

if depth == 0
    centers_xy = [];
    return;
end
nCols = size(binaryMask, 2);

[ygrid, xgrid] = ndgrid(1, -nCols:maxRadius*2:nCols);
centers_xy = [reshape(xgrid, [], 1), reshape(ygrid, [], 1)];

if depth > 1
    newCenters_xy = generateHexLattice(binaryMask, maxRadius, depth - 1);
    yOffset = maxRadius * sqrt(3); % this is fromm basic trigonometry of a 60deg angle

    % Duplicate the first row and adjust positions
    duplicateCenters_xy = newCenters_xy;
    duplicateCenters_xy(:, 1) = duplicateCenters_xy(:, 1) + maxRadius;
    duplicateCenters_xy(:, 2) = duplicateCenters_xy(:, 2) + yOffset;

    centers_xy = [centers_xy; duplicateCenters_xy];
end
end


function centers_xy = checkCircleCenterPoints(maxRadius, centers_xy, binaryMask)

dx = maxRadius;

% 1st: remove any circles exceeding the image's upper/lower extents
badcenters = dx+centers_xy(:, 1) > size(binaryMask, 2) | centers_xy(:,2)-dx < 1 |...
    dx+centers_xy(:,2) > size(binaryMask,1) | centers_xy(:,1)-dx < 1;

centers_xy(badcenters, :) = [];

% 2nd: index into binary image and remove any center points not on top of white  
badcenters2= zeros(size(centers_xy, 1), 1);
indices = sub2ind(size(binaryMask), round(centers_xy(:, 2)), round(centers_xy(:,1)));
badcenters2(:, 1) = binaryMask(indices); 
centers_xy(~logical(badcenters2), :) = [];  

% 3rd: create mask and test no elementwise multiplication with binary mask
P = size(centers_xy, 1);
area = zeros(P, 2);
for pp = 1:P
    mask = createMask(images.roi.Circle('Center',centers_xy(pp, 1:2), 'Radius', dx), binaryMask);
    area(pp, 1) = sum(mask(:));

    mask2 = mask & binaryMask;
    area(pp, 2) = sum(mask2(:));
end 
centers_xy(area(:,2) - area(:, 1) < -150, :) = [];

% 4th: quickly make sure you didn't yeet every option away 
assert(size(centers_xy, 1)>2)

end

function cleanedList = cleanedUpList(dirtyList)
% remove Inf's and NaN's

cleanedList = dirtyList;
  if any(isnan(dirtyList))
        nanidx = isnan(dirtyList(:,1)) | isnan(dirtyList(:,2));
        cleanedList(nanidx, :) = [];
  end
  if any(isinf(cleanedList))
        infidx = isinf(dirtyList(:,1)) | isinf(dirtyList(:,2));
        cleanedList(infidx, :) = [];
  end

end
