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

% In automated scripts, I usually will use try/catch when I call this
% function, because it will fail reasonably quickly in the event that it is not
% needed becase there are no gaps in the mask:

% try finalMask = LinkUpEdgeDiscontinuities(imcomplement(edgeImage));
% catch ME 
%     if strcmp(ME.identifier, 'MATLAB:badsubscript') 
%         finalMask = imcomplement(edgeImage);
%     else; rethrow(ME); end
% end

% Michael Glendinning, 2023
%%

% Find the candidate endpoints
eImage = imclose(edgeImage, strel(ones(13)));
branches = bwmorph(eImage, 'branchpoints');
theskel = eImage & ~branches;
theskel = bwareafilt(theskel, [90, Inf]);
theskel = imclose(theskel, strel(ones(5)));
theskel = imfill(theskel, 'holes');
theskel = bwskel(theskel, 'MinBranchLength', 110);

% "candidate" endpoints
endPs = bwmorph(theskel, 'endpoints');

% label image
[labeledImage, ~] = bwlabel(theskel);

labelvals = unique(labeledImage(find(endPs)));
newMaxDistances = struct('dist', [], 'index', [], 'ids', []);
for   h = 1:length(labelvals)
    val = labelvals(h); %labelvals is all the labels in thelabeled image ie each blob ID
    idx1 = find(labeledImage == val); %find every pixel with the same blob ID
    idx2 = find(endPs); %find each endpoint pixel
    ids = intersect(idx1 ,idx2); %identify the intersection (ie endpoints with the right label)
    
    for p = 1:numel(ids) %loop through the subset of the endpoints ids
        %bwdistgeodesic is a function that allows one to calculate distance constrained by mask. thus i get the "arc lengths" of along the skeleton, not as the "crow flies"
        D = bwdistgeodesic(theskel, ind2sub(size(edgeImage, 1:2), ids(p))); % calculate the geodescic distance
        D(D==Inf) = 0; %remove any any finities
        [newMaxDistances(h).dist(1, p), newMaxDistances(h).index(1,p)] = max(D, [], [1,2], 'linear'); %save to your structure which point was the farthest away from each point
    end
    
    [~,ii] = maxk(newMaxDistances(h).dist,2); %get the two points farthest from each other, these are your endpoints maximally far apart
    newMaxDistances(h).ids = ids(ii); %save the indices to the struct
    
end

%extract the coordinates of your endpoints from the array
maxApart_endPoints = struct2cell(newMaxDistances);
endPoints = cell2mat(squeeze(maxApart_endPoints(3,:,:)));
[endP_Rows, endP_Cols] = ind2sub(size(edgeImage), endPoints);

% "definitive endpoints"
endPoints_rc = [endP_Rows, endP_Cols];

% list all labels
allLabels = labeledImage(endPoints(1:numel(endPoints)));

%initialize loop
finalMask = edgeImage;

%you want to have this open so that the ROI createMask function knows the
%dimensions of the image to mask
f1 = uifigure('Visible', 'off');
ax1 = axes(f1);
imshow(edgeImage, [], 'Parent', ax1);
drawnow;
f1.Visible = 'off';

% special case! "u-shape" with only 1 gap in 1 blob, w/ 2 endpoints
cc = bwconncomp(edgeImage);
if cc.NumObjects == 1 && size(endPoints_rc, 1)==2
    hLine = images.roi.Line(ax1, 'Position',fliplr(endPoints_rc), 'LineWidth', 7);
    mask = createMask(hLine, edgeImage);
    finalMask = finalMask | mask;
    close all force
    return
end
% otherwise, we need to loop through the points to choose which points to
% connect to which...
for k = 1 : size(endPoints_rc, 1)
    thisPoint = endPoints_rc(k, 1:2); %[Y, X]
    thisLabel = allLabels(k);
    
    IndexSameLabeledPoints = (allLabels == thisLabel);
    
    viableOptions = endPoints_rc; % we don't want to connect points on the same blob!
    viableOptions(IndexSameLabeledPoints, :) = [] ;
    
    thepoints = [thisPoint; viableOptions];
    
    thepointsDist = pdist(thepoints, 'euclidean');
    distanceMat = squareform(thepointsDist);
    distances = distanceMat(1, 2:end);
    
    [~, minIdx] = min(distances);
    connectingPoint = thepoints(minIdx+1, :);
    
    % Draw line from this endpoint to the other endpoint. we have to
    % fliplr(ie. flip left - to - right) betcause we have [y, x] not [x,y]
    hLine = images.roi.Line(ax1, 'Position',[fliplr(thisPoint); fliplr(connectingPoint)], 'LineWidth', 7);
    mask = createMask(hLine, edgeImage);
    
    finalMask = finalMask | mask;	% Burn line into output image.
    
end
%finalMask = bwareafilt(finalMask, 1);

close all force



















