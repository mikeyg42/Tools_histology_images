 function [MOVINGmask_pad, MOVINGgray_pad, IMGmask_pad, IMGgray_pad] = pad4Images(minDistEdge, MOVINGmask, MOVINGgray, IMGmask, IMGgray)

%these images are the perimeter outline of their respective binary blobs
BW1 = bwperim(MOVINGmask);  
BW2 = bwperim(IMGmask);

% get image sizes. It is assumed MOVINGmask is same height and width as
% MOVINGgray, and likewise for IMGmask/IMGgray
[r1, c1]= size(MOVINGmask, 1:2);
[r2, c2]= size(IMGmask, 1:2);

%these images are all zeros except 1 pixel wide border of 1's around perimeter
perim1 = true([r1, c1]);
perim2 = true([r2, c2]);
perim1(2:end-1, 2:end-1) = false;
perim2(2:end-1, 2:end-1) = false;

% this calculates the distance between the point in the blob
% perimeter (ie BW1 or BW2) that is closest to one of the images edges. 
[d, ~] = pdist2(fliplr(find(perim1)), fliplr(find(BW1)), 'euclidean', 'smallest', 1); 
[e, ~] = pdist2(fliplr(find(perim2)), fliplr(find(BW2)), 'euclidean', 'smallest', 1);

% preset var "minDistEdge" sets up the closest distance the closest point
% in the image can be to an edge. minimumPad tells us how much more to ADD
% to ensure that.
minimumPad1 = minDistEdge - min(d(:));
minimumPad2 = minDistEdge - min(e(:));

% don't want any negatives!
pAmt1 = max(minimumPad1, 0);
pAmt2 = max(minimumPad2, 0);

% preallocate
MOVINGmask_pad = zeros([r1, c1]+[pAmt1*2, pAmt1*2], 'logical');
MOVINGgray_pad = ones([r1, c1]+[pAmt1*2, pAmt1*2], 'double');
IMGmask_pad = zeros([r2, c2]+[pAmt2*2, pAmt2*2], 'logical');
IMGgray_pad = ones([r2, c2]+[pAmt2*2, pAmt2*2], 'double');

MOVINGmask_pad(pAmt1:r1+pAmt1-1, pAmt1:c1+pAmt1-1) = logical(MOVINGmask);
MOVINGgray_pad(pAmt1:r1+pAmt1-1, pAmt1:c1+pAmt1-1) = im2double(MOVINGgray); 

IMGmask_pad(pAmt2:r2+pAmt2-1, pAmt2:c2+pAmt2-1) = logical(IMGmask);
IMGgray_pad(pAmt2:r2+pAmt2-1, pAmt2:c2+pAmt2-1) = im2double(IMGgray); 

 end