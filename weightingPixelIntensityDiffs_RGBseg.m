function outOfBoundsMask = weightingPixelIntensityDiffs_RGBseg(rgbIMG)
% This function makes great use of the efficient segmentation
% that is possible with "graydiffweight" function in the context of a
% homogenous background. 

% This function twice compares the intensity values of every 
% pixel to some reference value. Once using gray intensitys, once merging colors
% Using these weighted differences as the new intensity values for each pixel 
% makes seperating the fairly homogenous white flat background easy.  

% With images of approximate size = 5000x5000x3, this function runs in ~7.5 seconds 

% Michael Glendinning 2023

if ndims(rgbIMG)~=3
    disp('MUST BE RGB');
    return
end
 
%% step 1 - loop through R,G,B channels. Each pixel in each of 3xgrayscale ims is
%% assigned a weight depending on closeness in intensity to refVak in im(2,2,:). 
% then we log2-transform each image and then power transform im.^2 before
% recombining into an RGB image

rgbImage = zeros(size(rgbIMG, 1:3));
for ii = 1:3
    chl = rgbIMG(:,:,ii);
    W = graydiffweight(medfilt2(chl, [5,5]),2,2); %median filter is edge-aware! sometimes slow, but here we have time for it
    I = log2(W+eps); %add eps to avoid log(0) issue
    rgbImage(:,:,ii) = (1-rescale(I.^2, 0, 1)).^4;
end

M = size(rgbImage,1); N = size(rgbImage,2);
Q = M*N;

rgbStack = reshape(rgbImage, Q, 3);

%% step 2 - color segmentation using euclidean distances between pixel tristimulus vectors
% I arbitrarily am taking a region in the  center of the image that is
% sz = 60x60x3, which I then average for each channel and use that 1x3 vector as my
% benchmark against which I will evaluate the vector distances of every
% other pxl.

randomROI = rgbImage(round(M/2)-30:1:round(M/2)+30, round(N/2)-30:1:round(N/2)+30, 1:3);
avgColor = permute(mean(randomROI, [1 2]), [1, 3, 2]);

% Compute the Euclidean distance between all rows of image and our avgColor_vec. 
Dist_color = sqrt(sum((rgbStack - avgColor(ones(Q,1),:)).^2, 2));

%% step 3 - use otsu thresholding to quickly seperate dark from light pixels. 

%segment this distance plot using otsu thresholding
[counts,~] = histcounts(I,(0:0.05:1));
T_otsu = otsuthresh(counts);

% Initialize I (stacked form)
I = zeros(Q, 1, 'logical'); 

% apply thresh
I(Dist_color <= T_otsu ) = true; 
I(Dist_color > T_otsu ) = false;

% Reshape "I" back into M-by-N array.
I = reshape(I, M, N); 

%% step 4 -  morphologically "complete" the region mask 
%pick out the largest blob, and then fill any holes in the mask.
a = bwareafilt(I, 2);
outOfBoundsMask = imfill(a, 'holes');