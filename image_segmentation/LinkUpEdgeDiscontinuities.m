function finalMask = LinkUpEdgeDiscontinuities(edgeImage)
% finalMask = LinkUpBrokenEdges(edgeImage) INPUT MUST BE BW
% this function should connect any gaps in a black and white edge mask.
% WHITE = FOREGROUND
% black = background

% Procedure: This function first defines a smooth skeleton. Then, candidate endpoints
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
%%

testim = imread('/Users/jglendin/Dropbox - Michael/Dropbox/processedMSimages/MS_B_EGFL7_adjustedRGBImage_v1.tiff');
testim = imresize(testim, 0.25, {@oscResampling, 4});
imgPre = rgb2gray(preprocessRawRGBims(testim, 0));
edgeImage = edge(imgPre,'sobel');
edgeImage = bwareafilt(imfill(imclose(edgeImage, strel(ones(3))), 'holes'), 80);
% figure;
% imshow(edgeImage);
[boundaries, ~, N, ~] = bwboundaries(edgeImage); %traces all boundaries, exterior and holes
outerBoundaries = table(boundaries(1:N), cellfun(@(x) size(x, 1), boundaries(1:N), 'UniformOutput', false),...
    'VariableNames', {'Boundaries', 'nVerts'});

largeEnough = cell2mat(outerBoundaries.nVerts)>150;
realBoundaries = outerBoundaries(largeEnough,:) ;

cleanEdges = false(size(edgeImage, 1:2));
[sortedBoundys, ~] = sortrows(table2array(realBoundaries), 2, 'descend');
row1 = sortedBoundys{1};
cleanEdges(sub2ind(size(cleanEdges, 1:2), row1(1), row1(2))) = true;
k = convhull(row1);

circularity = zeros(size(realBoundaries,1),1, 'double');
meanDist = zeros(size(realBoundaries,1),1, 'double');
area = zeros(size(realBoundaries,1),1, 'double');
minDist = zeros(size(realBoundaries,1),1, 'double');

area(1,1) = polyarea(row1(k, 1), row1(k, 2));
perimeter = sum(sqrt(sum(diff(row1(k, :), 1, 1).^2, 2)));
circularity(1,1) = 4 * pi * area(1,1) / perimeter^2;

for curB = 2:size(realBoundaries,1)
    rowNext = sortedBoundys{curB,1};
    
    distances = pdist2(row1, rowNext, 'euclidean');   % compute the pairwise distances between X and Y
   % meanDist(curB, 1) = mean(mink(min(distances,[], 1), 80));
    minDist(curB, 1) = mean(mink(min(distances, [], 1), 4));
    k = convhull(rowNext);
    
    % Compute the area and perimeter of the convex hull, use to ]ntyd Compute the circularity
    area(curB, 1) = polyarea(rowNext(k, 1), rowNext(k, 2));
    perimeter = sum(sqrt(sum(diff(rowNext(k, :), 1, 1).^2, 2)));
    circularity(curB,1) = 4 * pi * area(curB,1) / perimeter^2;
    
    cleanEdges(sub2ind(size(cleanEdges, 1:2), rowNext(1), rowNext(2))) = true;
    row1 = [row1; rowNext];
end
data = [cell2mat(sortedBoundys(:,2)), minDist, meanDist, circularity, area];

sz = size(cleanEdges, 1:2);

cleanEdges1 = false(sz);
cleanEdges2 = cleanEdges1;
%cleanEdges3 = cleanEdges1;
cleanEdges4 = cleanEdges1;
cleanEdges5 = cleanEdges1;

[~, id1] = sortrows(data, -1); % nVerts
[~, id2] = sortrows(data, -2); % meanDist
%[~, id3] = sortrows(data, -3); % meanDist
[~, id4] = sortrows(data, 4 ); % Circularity
[~, id5] = sortrows(data, -5); % Area

for k = 1:round(size(realBoundaries, 1)/1.5)
    nextRow1 = sortedBoundys{id1(k), 1};
    nextRow2 = sortedBoundys{id2(k), 1};
   % nextRow3 = sortedBoundys{id3(k), 1};
    nextRow4 = sortedBoundys{id4(k), 1};
    nextRow5 = sortedBoundys{id5(k), 1};
    cleanEdges1(sub2ind(sz, nextRow1(:, 1), nextRow1(:, 2))) = true;
    cleanEdges2(sub2ind(sz, nextRow2(: ,1), nextRow2(:, 2))) = true;
    %cleanEdges3(sub2ind(sz, nextRow3(:, 1), nextRow3(:, 2))) = true;
    cleanEdges4(sub2ind(sz, nextRow4(: ,1), nextRow4(: ,2))) = true;
    cleanEdges5(sub2ind(sz, nextRow5(: ,1), nextRow5(: ,2))) = true;
    drawnow;
end
combIm = cat(3, cleanEdges1, cleanEdges2 , cleanEdges4, cleanEdges5);
combIm = sum(im2double(combIm), 3)./4;
myboundary =  combIm > 0.7;
% figure;
% imshow(myboundary);

[i, j] = find(myboundary);
Bcoords = boundary([i, j], 0.65);
y = i(Bcoords);
x = j(Bcoords);
%hold on
%plot(x,y,'r-o');
%hold off
mask = poly2mask(x, y, size(cleanEdges, 1), size(cleanEdges, 2));
mask = mask & ~imerode(mask, strel(ones(3)));

% figure;
THEboundary = imreconstruct(edgeImage&mask, edgeImage);
THEboundary =imclose(imfill(THEboundary, 8, 'holes'), strel('diamond',3));
% imshow(THEboundary);

theskel = bwskel(THEboundary, 'MinBranchLength', 40);
branches = imdilate(bwmorph(theskel, 'branchpoints'), strel(ones(3)));
theskel = theskel & ~branches;

endPs = bwmorph(theskel, 'endpoints');
[y,x] = find(endPs);

D = squareform(pdist([x, y]));
D(find(eye(size(D)))) = [];
D = reshape(D, [], size(x,1)-1);
[a,b] = min(D, [],2);

[labeledImage, ~] = bwlabel(theskel);
labelvals = unique(labeledImage(sub2ind(sz, y, x)));

endPointArray = cell(numel(labelvals), 1);
for   nLoops = 1:length(labelvals)
    ePoints = zeros(2,2);
    val = labelvals(nLoops); %labelvals is all the labels in thelabeled image ie each blob ID
    idx1 = find(labeledImage == val); %find every pixel with the same blob ID
    idx2 = find(endPs); %find each endpoint pixel
    ids = intersect(idx1 ,idx2); %identify the intersection (ie endpoints with the right label)
    
    if numel(ids)>2
        disp('theres a problem');
    else
        [ePoints(1, 1), ePoints(1, 2)] = ind2sub(sz, ids(1));
        [ePoints(2, 1), ePoints(2, 2)] = ind2sub(sz, ids(2));
        endPointArray{nLoops, 1} = ePoints;
    end
end


% create endPoint map
endPoints = containers.Map(labelvals,endPointArray,'UniformValues',false);

% <> unfinished </>





%initialize loop
finalMask = THEboundary;

%you want to have this open so that the ROI createMask function knows the
%dimensions of the image to mask
f1 = uifigure('Visible', 'on');
ax1 = axes(f1, 'Visible', 'on');

imshow(finalMask, 'Parent', ax1);
drawnow;

% This is in the rare event you just have 1 blob with a sole discontoinuity (a "U" blob)
if length(labelvals)== 1
    hLine = images.roi.Line(ax1, 'Position',fliplr(ePoints), 'LineWidth', 7);
    mask = createMask(hLine, edgeImage);
    finalMask = mask | finalMask;
    close all force
    return
end


% first lets draw in the really tiny lines. 
% then with the remaining points lets either draw all the lines or all the lines interior
% or the first closest lines, 
% then we can delete the lines that are bad until we are happy 



%=============
finalMask = cleanMask(finalMask); %this function I wrote adaptively size filters blobs

% DO THE IMFILL CHECK BEFORE GIVING UP!


close all force



















