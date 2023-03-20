function [fixedMat_new,movingMat_new, myNEWimages] = preReg_padding_center(fixedMat, movingMat, myOLDimages)
% [fixedMat_new,movingMat_new, myNEWimages] = preReg_padding_center(fixedMat, movingMat, myOLDimages)

% note that: myOLDimages = {MOVINGmask2, MOVINGgray2, ref_IMGmask,
% ref_IMGgray}; and myNEWimages has the same order in its cell array

% This is not desgined for RGB images or anything w/ third dimension. just use
% imsplit in that event and feed channels though
%
%fixed and moving mat must be n x 2 [x,  y] coordinates. This was written
% with the intention of centering my quadrilaterials' corner points + centroid 
% atop one another. But will actually work for any size x-by-2 matrix of
% coordinates from the same images you're registering 
%
% all fixed images must be the same size as other fixed ims!
% all moving images must be the same size as the other moving ims!! 
% any logical masks

%% first we adjust the matrices so that the arithmetic mean of all coordinates 
%% in one matrix is on top of the arithmetic mean of the other set of coordinates
%% but by only moving into the positive direction

diff = fixedMat - movingMat;
avg = sum(diff, 1)./size(diff, 1);

positive_transl = avg; 

% PROBLEM: we could have just translated movingMat this amount [dx, dy] and successfully 
% aligned these pointclouds. But if dx or dy < 0, our images will become
% clipped, which, sure, is no big deal for a matrix, but if images get trucated,...
% its BAD!

% SOLUTION: split up "+" and "-" displacements, and in event displacement<0, apply 
% the inverse, "+" translation to other image. 
%  ... w/all displacement going in a "+" direction , nothing is ever clipped
% example: 
% avgDiff = [+10, -8]. in the X direction, movingIm displaced [0, +10].
% for Y displacement, we shift the FixedIm [0, +8] (and NOT movingIM [0,
% -8]). 

%These shifts are not translations, but instead precise padding of
% image borders.

negative_transl = positive_transl.*-1;
positive_transl(positive_transl<0)=0;
negative_transl(negative_transl<0)=0;

%now we apply translation to the matrices:
movingMat_new = movingMat + repmat(positive_transl, size(movingMat,1), 1)+1; %we need to add +1 because we are going to add the same to the images shortly
fixedMat_new = fixedMat + repmat(negative_transl, size(fixedMat,1), 1)+1;%we need to add +1 because we are going to add the same to the images shortly

%% turning now to our images...
% although we above could maintain subpixel accuracy w coordinate
% matrices, padarray needs integers. padarray also needs [row, col], not [x, y]
negative_transl = round(fliplr(negative_transl));
positive_transl = round(fliplr(positive_transl));

% "translate" images using padarray to align centroids as you did mats above
movingIms2_1 = padarray(myOLDimages{1}, positive_transl, 0, 'pre'); % mask
movingIms2_2 = padarray(myOLDimages{2}, positive_transl, 1, 'pre');
fixedIms2_1 = padarray(myOLDimages{3}, negative_transl, 0, 'pre'); %mask
fixedIms2_2 = padarray(myOLDimages{4}, negative_transl, 1, 'pre');

fixedIms2 = cat(3, fixedIms2_1, fixedIms2_2);
movingIms2 = cat(3, movingIms2_1, movingIms2_2);

% registration sometimes demands equal sized moving and fixed images... 
% so: make a square canvas w/ side length equal to the longest side of all the images
[r1, c1] = size(fixedIms2, 1:2);
[r2, c2] = size(movingIms2, 1:2);
mx = max([r1, r2, c1, c2]);

newIm1 = cat(3, zeros(mx+2, mx+2, 2, 'double'), ones(mx+2, mx+2, 2, 'double'));
newIm2 = cat(3, zeros(mx+2, mx+2, 2 ,'double'), ones(mx+2, mx+2, 2, 'double'));
%by adding this plus 2, we in effect give an extra 1 pixel of padding, which
%should prevent any issues with later on calling imclearborders.

%adjusting the matrices accordingly 
%bc we added +2 to the image, thats a shift of +1 away from origin, and +1 extent increase on the periphery. \
% The matrix only cares about the shift, so we ignore the later +1.
movingMat_new = movingMat_new + repmat(1, size(movingMat_new,1), 2); 
fixedMat_new = fixedMat_new + repmat(1, size(fixedMat_new,1), 2);

% assign images to the same spot to stay in alignedment but one away from the edge so no edge issues should arise. 
newIm1(2:r1+1, 2:c1+1, 1:2) = fixedIms2;
newIm2(2:r2+1, 2:c2+1, 1:2) = movingIms2;

% split channels back and make into 1 cell array for final step
[fixedIm1mask, fixedIm2p] = imsplit(newIm1);
[movingIm1mask, movingIm2p] = imsplit(newIm2);
myNEWimages = {movingIm1mask, movingIm2p,fixedIm1mask, fixedIm2p};

% final step is to reformat and fix up the masks!!
mask1 = imclearborder(logical(myNEWimages{1}));
mask2 = imclearborder(logical(myNEWimages{3}));

mask2 = imfill(mask2, 'holes');
mask1 = imfill(mask1, 'holes');

myNEWimages{1} = mask1;
myNEWimages{3} = mask2;

end



