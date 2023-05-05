function finalMask = LinkUpEdgeDiscontinuities(inputImage)
% finalMask = LinkUpBrokenEdges(BWedgeImage)
% OR
% finalMask = LinkUpBrokenEdges(BWedgeImage, RGBimage, parentFunctionCall)
% Varargin{1} is RGB im. Varargin{2} feeds the function info on the parent fucntion
% calling it.

% This function should connect any gaps in a black and white edge mask.
% WHITE = FOREGROUND
% black = background

% We give the option to incorporate RGBimage data and or information about the parent
% function calling runLevelSets.

% Procedure proper: This function first defines a smooth skeleton. Then, candidate endpoints
% are selected using bwmorph(I, 'endpoints'). Subsequently, we use BWLABEL
% to label each contiguous blob in the image. Because there are usually
% many more endPoints on a blob than we want (which is here set to be
% 2/blob), we have nested for-loops that go through each each blob, pull
% the endpoint candidates they contain, and then loop through those
% endpoints, calculating distances between each via BWDISTGEODESIC in order to find the
% 2 points maximally far apart as constratined by the mask of the blob.
% This leaves us with 2 end points per blob. Then the function goes through
% each of these definitive end points and connects it to the end point
% closest to it that is on a different blob. A straightline ROI draw,
% converted to a mask, which is added to edgeImage.

% Michael Glendinning, 2023
% ======================================================================== %
% inputImage = inputImage.*imcomplement(compositeMask);
% close all
% 
% if ~islogical(inputImage) && ~(all(unique(inputImage(:))==0|unique(inputImage(:))==1) && isfloat(inputImage))
%     error('format of inputImage must be BINARY. either logical, or a floataing precision image with ONLY 0s and 1s. ')
% end
% 
% if nargin>2
%     parentFunction = varargin{2};
% else
%     parentFunction = 'none';
% end
% 
% if contains(parentFunction, 'Level', 'IgnoreCase',true)
%     inputImage_thin = bwmorph(inputImage,'thin');
%     edgeImage = bwmorph(bwmorph(inputImage_thin, 'spur'), 'spur');
%     edgeImageSK = bwskel(edgeImage);
%     if checkImfill(edgeImageSK)
%         finalMask = edgeImageSK;
%         disp('there were no discontinuities in the LevelSet result. Returning to runLevelSet');
%         %return
%     end
% else
%     edgeImage = inputImage;
%     edgeImageSK = bwskel(edgeImage);
% end
% 
% rgbFlag = 0;
% if nargin>1 && ndims(varargin{1})==3 %then we have an RGB image
%     rawRGBim = im2single(varargin{1});
%     RGBim = imresize(rawRGBim, size(inputImage, 1:2));
%     rgbFlag = 1;
% end
% 
% 
% if ~checkQuadrants(edgeImageSK) && rgbFlag == 0 %if you have the rgbFlag, then we can try another segmentation/edge detection
%     disp('you are missing a quadrant? unsalvagable');
%     return;
% elseif ~checkQuadrants(edgeImageSK) && rgbFlag == 1
%     [~, BWedge] = fullColorSegmentation(RGBim);
% 
%     if checkImfill(BWedge)
%         finalMask = BWedge;
%         return;
%     elseif checkImfill(BWedge|edgeImage)
%         finalMask = BWedge|edgeImage;
%         return;
%     elseif checkQuadrants(BWedge|edgeImage)
%         edgeImage = BWedge|edgeImage;
%         edgeImageSK = bwskel(imclose(edgeImage, strel(ones(7))));
%     elseif ~checkQuadrants(BWedge|edgeImage)
%         error('unable to draw any continuous edge');
%     end
% end
%%

testim =  inputImage;
testim = imresize(testim, 0.25, {@oscResampling, 4});
labim = rgb2lab(testim);
%grayIm = imadjust(rescale(labim(:,:,1), 0,1));
grayIm = histeqfloat(rescale(labim(:,:,1), 0, 1));


% find the local variance and then threshold it
localVarianceMap = stdfilt(grayIm, true(3));
T = multithresh(localVarianceMap, 5);
% use morphological reconstruction
[dx, dy, dxy, dxx, dyy] = derivative7(imdiffusefilt(grayIm), 'x', 'y', 'xy','xx', 'yy');
SGDG = (dxx.*dx.^2+dyy.*dy.^2+2.*dxy.*dx.*dy)./(dx.^2+dy.^2);

edgeMap = imreconstruct(edge(SGDG, 'log'), imbinarize(localVarianceMap , T(3)));

fractDeriv1 = imFractionDeriv(grayIm, 0.5);
fDeriv1 = locallapfilt(im2single(fractDeriv1), 0.5, 1.5);


fractDeric2 = imFractionDeriv(grayIm, 1.5);
fDeriv2 = locallapfilt(im2single(fractDeric2), 0.8, 1.5);


% taking the first 500 largest blobs should be PLENTY
edgeIm = imclose(edgeMap, strel(ones(7)));
bigEdgeIm = bwareafilt(edgeIm, 50);
smallEdgeIm = bwareafilt(edgeIm, 500);
smallEdgeIm(bigEdgeIm == 1) = 0;
bigBlobs = bwperim(imclose(bigEdgeIm, strel(ones(15))));
[y,x] = find(bigBlobs);
bds = boundary([x,y], 0.5);
perim = bwperim(poly2mask(x(bds), y(bds), sz(1), sz(2)));
perim2 = imdilate(perim, strel('octagon', 132));
%% Start of Algorithm proper

sz = size(edgeImageSK, 1:2);

% calculte data for EACH blob in the bw image's boundary using BWBOUNDARIES, make table

[~, blobLabels] = bwboundaries(edgeImageSK);


c = regionprops('table', allEdgeIm, 'Centroid');
centersD = squareform(pdist(c.Centroid, 'fasteuclidean'));
centersD = sort(centersD, 1);
centersD()

scoreMatrix = zeros(max(blobLabels(:)));


       


outerBoundaries = [arrayfun(@(x) size(boundaries{x,1},1), 1:1:size(boundaries, 1)); 1:1:size(boundaries, 1)]' ;

% only keep blobs longer than 200 pixels
largeEnough = cell2mat(outerBoundaries.nVerts)>200;
realBoundaries = outerBoundaries(largeEnough,:) ;

% initallize some structures for data in upcoming loop that assesses each blob
area = zeros(size(realBoundaries,1),1, 'double');
circularity = zeros(size(realBoundaries,1),1, 'double');
minDist = zeros(size(realBoundaries,1),1, 'double');

% <p> note that we had to start extracting data for row 1 outside the loop because minDist
% is recursively defined...for the 2nd blob, its "distance" is measured against only
% blob #1. Then blob1 and blob2 "merge" and blob3 is measured against them both...
% and so on. This counteracts the disanvantage of small but potentially critical blobs </p>

[sortedBoundys, ~] = sortrows(table2array(realBoundaries), 2, 'descend');
row1 = sortedBoundys{1};
k = convhull(row1);
area(1,1) = polyarea(row1(k, 1), row1(k, 2));
perimeter = sum(sqrt(sum(diff(row1(k, :), 1, 1).^2, 2)));
circularity(1,1) = 4 * pi * area(1,1) / perimeter^2;

for curB = 2:size(realBoundaries,1)
    rowNext = sortedBoundys{curB,1};

    distances = pdist2(row1, rowNext, 'euclidean');   % compute the pairwise distances between X and Y
    % meanDist(curB, 1) = mean(mink(min(distances,[], 1), 80));
    minDist(curB, 1) = mean(mink(min(distances, [], 1), 4));
    k = convhull(rowNext);

    % Compute the area and perimeter of the convex hull, use to Compute the circularity
    area(curB, 1) = polyarea(rowNext(k, 1), rowNext(k, 2));
    perimeter = sum(sqrt(sum(diff(rowNext(k, :), 1, 1).^2, 2)));
    circularity(curB,1) = 4 * pi * area(curB,1) / perimeter^2;

    row1 = [row1; rowNext];
end
% Extract blob data and convert to matrix
data = [cell2mat(sortedBoundys(:,2)), minDist, circularity, area];

[~, id1] = sortrows(data, -1); % nVerts
[~, id2] = sortrows(data, -2); % minDist
[~, id3] = sortrows(data, 3); % Circularity
[~, id4] = sortrows(data, -4); % Area

id1(round(size(realBoundaries, 1)/1.5):end) = [];
id2(round(size(realBoundaries, 1)/1.5):end) = [];
id3(round(size(realBoundaries, 1)/1.5):end) = [];
id4(round(size(realBoundaries, 1)/1.5):end) = [];

[counts, edges] = histcounts([id1; id2; id3; id4], 1:1:size(realBoundaries, 1)+1);
counts = counts./4;

combinedImg = zeros(sz, 'double');
for k = 1:size(realBoundaries, 1)

    combinedImg(realBoundaries(k ,1).Boundaries) = counts(1, k);
end


for k = 1:round(size(realBoundaries, 1)/1.5)
    nextRow1 = sortedBoundys{id1(k), 1};
    nextRow2 = sortedBoundys{id2(k), 1};
    nextRow3 = sortedBoundys{id3(k), 1};
    nextRow4 = sortedBoundys{id4(k), 1};

    cleanEdges1(sub2ind(sz, nextRow1(:, 1), nextRow1(:, 2))) = 1;
    cleanEdges2(sub2ind(sz, nextRow2(: ,1), nextRow2(:, 2))) = 1;
    cleanEdges3(sub2ind(sz, nextRow3(: ,1), nextRow3(: ,2))) = 1;
    cleanEdges4(sub2ind(sz, nextRow4(: ,1), nextRow4(: ,2))) = 1;
end

% Combine binary images and convert to grayscale image
combIm = im2double((cleanEdges1 + cleanEdges2 + cleanEdges3 + cleanEdges4) ./ 4);

% Assess the histogram of nonzero values
[counts, histedges] = histcounts(combIm(combIm~=0));
cumDist = cumsum(counts)/sum(counts(:));

% threshold distribution such that there is no barren image quadrant and not less than
% a third of the non-zero pixels
test = false; % initialize
t = 3; % begin at 1/3
while ~test && t > 1
    T = round(1/t, 4); % e.g. as t decreases, T -> 0.33, 0.40, 0.50, 0.67
    cumD = cumDist<T;
    % each val is a cutoff thresh and the % of im kept [1, 1, ...,1, 1, 0, 0,...,0 ]
    edgeIndex = int8(sum(cumD(:))-1);
    if edgeIndex>0 && isinteger(edgeIndex)
        TT = histedges(edgeIndex);
    else %not sure why MATLAB made me do this but its superfluous, edgeIndex will always be a positive integer
        edgeIndex = round(abs(edgeIndex));
        TT = histedges(edgeIndex);
    end
    myboundary = imbinarize(combIm, TT) ;     % Threshold at that level
    test = checkQuadrants(myboundary); % see if basic quadrant test is passed
    if ~test
        t = t-0.35;  % ... if it isnt passed, lower t, thereby raising T
    end
end

% preset structuring elements
se5 = strel(ones(5)); se3 = strel(ones(3));

% Use imfill cautiously. then after exclude again using size
if checkImfill(myboundary)
    edgeI = imfill(myboundary, 'holes');
else
    edgeI  = imclose(myboundary, se5);
end
edgeI = bwareaopen(edgeI, 200);

% Intention here is to transform the original image such that the inside is flat black suych that it will "punch
% out the reamining of white in the middle of the image

    grayTestim = locallapfilt(im2single(grayIm), 0.6, 3); %heavy, edge-aware filtering!
    % preset loop
    nLevels = 10;
    T = multithresh(grayTestim, nLevels);
    areaPostClear = zeros(nLevels, 1);
    for lvl = 3:nLevels
        invThresh = imcomplement(imbinarize(grayTestim, T(lvl)));
        testim = imclearborder(invThresh);
        areaPostClear(lvl, 1) = sum(testim(:));
        if areaPostClear(lvl,1)*15 < areaPostClear(lvl-1, 1)
            break
        end
    end

% The image we want will be where the areaPostClear value is maximized. 
[~, Tidx] = max(areaPostClear);
% By taking only largest blob in inverse image, we get just the background around the black blob
outsideMask = bwareafilt(imbinarize(grayTestim, T(Tidx)), 1);
 
edgeMap = imdilate(edgeI, se5) & ~imcomplement(outsideMask) ;
distMap = bwdist(imclose(edgeMap, se5));
thickEdge = imcomplement(bwmorph(distMap, 'skel'));
avgEdge1 = im2double(edgeMap)+im2double(thickEdge);

%% incorporate fractional derivatives!
fractDeriv1 = imFractionDeriv(grayIm, 0.5);
fDeriv1 = locallapfilt(im2single(fractDeriv1), 0.5, 1.5);
avgIm1 = rescale(imreconstruct(im2single(avgEdge1), imnlmfilt(fDeriv1)),0,1);

fractDeric2 = imFractionDeriv(grayIm, 1.5);
fDeriv2 = locallapfilt(im2single(fractDeric2), 0.8, 0.2);
fDeriv2 = histeqfloat(fDeriv2).^4;

avgIm2 = imclose(imcomplement(avgIm1), se5) + imclose(fDeriv2, se5);

derivEdge = bwperim(imbinarize(rescale(avgIm2,0,1)));
avgEdge2 = derivEdge & bwareaopen(imbinarize(avgIm2), 500);
THEedge = im2double(bwareaopen(avgEdge2, 400));

if checkImfill(THEedge)
THEmask = imfill(THEedge, 'holes');
finalMask = bwareafilt(THEedge,1);
else 
    % maybe its just a tiny gap that needs filling - try to close image
    THEedge = imclose(THEedge, se5);
    
    if checkImfill(THEedge)
    THEmask = imfill(THEedge, 'holes'); 
    finalMask = bwareafilt(THEedge,1);

    else 
        finalMask = false(sz);
    end
end
end

%if sum(THEmask(:)) == 0
% % Re-query new skeleton for endpoints
% endPs = bwmorph(THEedge, 'endpoints');
% [y,x] = find(endPs);
% endPointCoords = [x,y];
% endPointIndx = find(endPs);
% 
% % Use PDIST to measure distances bewteen points.
% D = squareform(pdist(endPointCoords));
% D(eye(size(D))==1) = 9999;
% [dists, distInd] = min(D, [],2); % use @min to determine the closest point to each endP
% 
% % When the nearest neighbor to an endpoint is within 7pixels away, we can be fairly
% % confident that the line connecting them will help us
% veryShort = dists < 11;
% shortIndex = [find(veryShort), distInd(veryShort)]; % left column has the owners fo the cols with min <11
% % right col has the index of what the left was closest TO
% % loop through these short points and draw the tiny lines between points
% figure; imshow(THEedge);
% 
% for sh = 1:size(shortIndex,1)
%     pt1 = shortIndex(sh, 1); % pt1/2 gives the INDEX
%     pt2 = shortIndex(sh, 2);
%     point1 = [x(pt1), y(pt1)]; % point1/2 is reformatted as xy coords
%     point2 = [x(pt2), y(pt2)];
%     roi = drawline('Position', [point1; point2]);
%     roiMask = createMask(roi);
%     finalSkel = finalSkel | imdilate(roiMask, se3); % dilate bc the lines otherwise will be 1 pixel thick
% end
% close all
% finalSkel = bwskel(finalSkel);
% 
% 
% % Now we must connect the points that are farther apart than 10 pixels...
% % to do we first just identify our "true" end points
% notdonepoints = setdiff(1:1: numel(dists), unique(reshape(shortIndex, [], 1), 'sorted'));
% lastEndPs = zeros(sz, 'logical');
% lastEndPs(endPointIndx(notdonepoints)) = 1;
% idx1 = find(lastEndPs); %find each endpoint pixel
% 
% % we don't want to draw lines between points that share a BWLABEL. 
% [labeledImage, ~] = bwlabel(finalSkel);
% labelvals = unique(labeledImage(idx1)); % idx1 restricts to the labels containing an endpoint
% 
% if numel(labelvals) == 2
%     % connect your U-shape
% 
% else
% % Pre-allocate cell array for results of foor-loop
% endPointArray = cell(numel(labelvals), 1);
% 
% for nLoops = 1:length(labelvals) %loop through each label from BWLABEL
%     ePoints = zeros(2,2); % Preallocate/clear matrix for endpoint coords
%     val = labelvals(nLoops); %labelvals is all the labels in the labeled image held by at least 1 endpoint identified
%     idx2 = find(labeledImage == val); %find every pixel with the same blob ID
% 
%     %identify the intersection (ie the endpoints sharing the particular loop iteration's label val
%     ids = intersect(idx1 ,idx2);
% 
%     [yi, xi] = ind2sub(sz, ids);
%     coords = [xi, yi];
%     nMatches = numel(ids);
%     if nMatches>2 %CASE 1: (unlikely) 1 blob has multiple points we've tagged as maybe endpoints
% 
%         D = zeros(nMatches, nMatches);
%         for p= 2:nMatches
%             for q = 1: nMatches
%                 if p>q
%                     D(p,q) = DistFun(finalSkel, coords(p, 1:2), coords(q, 1:2));
%                 end
%             end
%         end
%         [f,g] = find(max(D [], "all"));
% 
%         [ePoints(1, 1), ePoints(1, 2)] = ind2sub(sz, ids(f)); % After exhaustive search of potential points of a blob, bestPair is furtherest 2 apart
%         [ePoints(2, 1), ePoints(2, 2)] = ind2sub(sz, ids(g));
% 
%     elseif numel(ids) == 1 % I don't entirely know yet why this might happen but it can happen clearly
% 
%         [labelIm,~] = bwlabel(finalSkel);
%         thisval = labelIm(ids);
% 
%         [endpIx, ~] = find(endPointIndx == ids);
%         thisCoord = endPointCoords(endpIx, 1:2); %xy
% 
%         othervals = false(sz);
%         othervals(labelIm ~= thisval & labelIm~=0) = true;
%         [i, j] = find(othervals);
%         D = pdist2([j,i], thisCoord, "euclidean");
%         [~, iix] = sort(D, 'ascend');
%         ix = iix(1);
%         ePoints(1, 1) = thisCoord(1);
%         ePoints(1, 2) = thisCoord(2);
%         ePoints(2, 1) = j(ix);
%         ePoints(2, 2) = i(ix);
% 
%     else % CASE 2: (most likely) 1 blob has only 2 points identified, which we assume are perfect, exemplar endPoints.
%         [ePoints(1, 1), ePoints(1, 2)] = ind2sub(sz, ids(1));
%         [ePoints(2, 1), ePoints(2, 2)] = ind2sub(sz, ids(2));
%     end % we should NEVER see the case of 0 points because we defined our set using the intersection of endpoint and shared label defining loop
% 
%     endPointArray{nLoops, 1} = ePoints; % output of cases 1 and 2 are efficiently stored for later
% end
% 
% % use nearly filled skeleton image to approximate image with a polygon.
% [yf, xf] = find(THEedge);
% pgon = polyshape(xf, yf);
% 
% % We leverage the polyshape classes' function "polybuffer" to draw a buffer around
% % image just a few pixels wide, which we will use to mask away extraneous solutions
% figure; imshow(THEedge);
% pgon2 = polybuffer(pgon, 5.5, 'JointType','square');
% hold on; plot(pgon2, 'FaceAlpha', 0.5);
% bufferMask = createMask(pgon2);
% exteriorChopper = imcomplement(bufferMask);
% mask2 = padarray(createMask(pgon), [5,5], 0,'both'); 
% interiorChopper = imresize(mask2, sz, {@oscResampling, 4});
% chopper = exteriorChopper | interiorChopper;
% %
% figure;
% imshow(finalSkel); hold on; visboundaries(chopper);
% blankIm =  false(sz);
% for L = 1:size(endPointArray, 1)
% 
%     pt1s = cell2mat(endPointArray(L)); %xy
%     pt2s = cell2mat(endPointArray(labelvals>L));  % this inequality cut the time of this loop in half by pop[ulating only half of the symmetrical triangular matrix
% 
%     % for each label value, we now have two points (pt1s andn pt2s). We will now draw EVERY
%     % line possible starting at these two spots, and going ALL the spots with DIFFERENT labels
%     for p = 1:size(pt2s, 1)
%         line1 = [pt1s(1, :); pt2s(p, :)];
%         % "line 1"
%         roi1 = drawline('Position', fliplr(line1));
%         r1 = createMask(roi1);
%         test1 = chopper.*r1;
%         if any(test1(:)~=0) 
%             line1 =  blankIm;
%         else
%             line1 = imdilate(r1, se3);
%         end
%         % "line 2"
%         line2 = [pt1s(2, :); pt2s(p, :)];
%         roi2 = drawline('Position', fliplr(line2));
%         r2 = createMask(roi2);
%         test2 = chopper.*r2;
%         if any(test2(:)~=0) 
%             line2 =  blankIm;
%         else
%             line2 = imdilate(r2, se3);
%         end
% 
%         % dilate both masks slightly and then  burn them ALL into the final skeleton
%         finalSkel = finalSkel | line1 | line2;
%     end
% end
% figure;
% imshowpair(finalSkel, finalSkel2);
% 
% finalEdge  = finalSkel.* imcomplement(bufferMask);
% 
% if sum(sum(bwperim(filledImage)))/sum(finalSkel(:)) > 1.8 || ~checkImFilled(finalEdge)
%     disp('final mask has a problem....');
%     if rgbFlag == 1
%         disp('returning edge boundary generated with a different segmentation algorithm using the provided RGB image')
%         [~, finalMask] = bwperim(fullColorSegmentation(RGBim));
%     end
% else
%     finalMask = imclose(finalEdge, strel('disk', 9, 8));
%     if ~checkImFilled(finalMask)
%         error('not able to remove the discontinuities')
%     elseif rgbFlag == 1
%        im = rgb2lab(RGBim);
%        avgIm2 = imadjust(rescale(im(:,:,1), 0, 1));
%        finalMask = activecontour(avgIm2, finalMask, 15);
%     end
% 
% end     
%end
%end

