function outOfBoundsMask = hysteresisThreshold_wMorph(grayImage, RGBimage)
% syntax =  outOfBoundsMask = hysteresisThreshold_wMorph(grayImage, RGBimage)

% my implementation of a hysteresis type of thresholding (also known as
% double threshholding). I've opted to use histogram (i.e. spectral)
% attributes to set up the lower threshold level (the "wide net"). And then
% I use IMEXTENDMAX, which calculates the H-maxima transform to locate
% definitive peaks, which I can then use to index into the wide-net segmentation to 
% remove extraneous object, which is accomplished with IMFILL and


%% step 1: evaluate the histogram and determine otsu's threshold level
GT = graythresh(grayImage); %gives the threshold level (via Otsu's method)

[n,bins] = imhist(RGBimage);
histminima = islocalmin(n,'MaxNumExtrema', 4); % should give all of the local minimums
minVals = bins(histminima); % what we want to do is locate that minimum nearest GT
[~,b]= min(abs(GT-minVals)); % find the closest local minimum to the gray thresh value

%% step 2: Cast "Wide-net" threshold
GenerousThreshold = minVals(b)*0.9; % just to be surely maximally inclusive we try going even more generous. too much and we risk it all going white
checkIt = imbinarize(grayImage, GenerousThreshold);
chk = sum(~checkIt(:));
if chk < 2000  %make sure you didn't overcompensate threshold and make it 100% white
    GenerousThreshold = minVals(b)*0.95;
    checkIt= imbinarize(grayImage, GenerousThreshold);
    chk = sum(~checkIt(:));
    if chk < 2000
        GenerousThreshold = minVals(b);
    end
end
checkIt= imbinarize(grayImage, GenerousThreshold);
chk = sum(~checkIt(:));

if sum(sum(imclearborder(imbinarize(grayImage, minVals(b))))) == chk
    if sum(sum(imclearborder(imbinarize(grayImage, GenerousThreshold))))<5000
        GenerousThreshold = minVals(b);
    end
end
% the above asks: "If, at the thresh. level of the histogram min, then
% there is part of mask touching border, then ok fine whatever - we
% just move on. BUT, we could get lucky and discover a border touching (i.e. after
% calling imclearborder its area is the same as chk). then, we can
% as if we've pushed GenerousThreshold TOO low another way. If that
% lower threshold DOES touch border significantly where the original histogram
% min level thresholding did not, that's maybe way too much.

%% Step 3: Complete Hysteresis with tighter net
% First we use strict thresh to define our seed areas (after a bit more
% cleaning of the images with clear borders and imadjust)
maxAreas = imextendedmax(imadjust(imclearborder(grayImage)), 0.5, 8);
maxArea_1 = bwareafilt(maxAreas, 1);
maxSeedIdx = find(maxArea_1);

% then we use imfill function to select the continuous areas that contain the seed areas within "generouser area")
generousBinary = imbinarize(grayImage, GenerousThreshold);
hysteresisImage = imfill(~generousBinary, maxSeedIdx, 8); % connectivity = 8
hysteresisImage = hysteresisImage & generousBinary;

%% postscript

% Just in case there are any little holes, we do a image close quickly
se = strel('disk', 9, 8);
outOfBoundsMask = imclose(hysteresisImage, se);

figure; imshow(outOfBoundsMask);

