function binaryMask = morphologySegment(imAdjRGB)
% imAdjRGB MUST be a double RGB image (3 channels only, no alpha chl). 
% It is assumed that background correction (I use the complement of the tophat of
% the complement image along with imadjust) has been completed already. Uses a whole lot
% of morphological functins to segment an RGB image into a binary foreground/background
% mask.
% Michael Glendinning, 2023

se1 = strel('diamond', 11);

compIm = ones(size(imAdjRGB,1:3));
compIm = compIm - imAdjRGB;
% 1-RGB is same as imcomplement(imAdjRGB), but prefer this bc maybe a bit faster, also imcomplement might rescale, which I dont want

compIm = imclose(compIm, se1); % morphological manipulation. because SE is 2d, it automatically treats each channel as a grayscale image
sumIm = sum(compIm,3);  %sum together all of the channels to create one image that should approximately span [0, 3] in intensity.

sumIm(sumIm>1) = 1; %clip the values above 1 as 1. 
sumIm(sumIm<0) = 0; %clip any negative values as 0's

% to suppress the brightest pixels, I create an im filtered with a local minimum filter convolution,
minfilt = @(x) min(x,[],1);
m = 7; n = 7; 
sumIm_padded = padarray(sumIm,[m, n],'replicate');
sImFilt = colfilt(sumIm_padded,[m n],'sliding', minfilt); %colfilt appears to be the fastest way to do this math!! but, I've been warned it consumes much memory!!
sImFilt = sImFilt(m+1:size(sImFilt,1)-m, m+1:size(sImFilt,2)-m); %remove padding

% its too DULL now... so, we can locate just the top 3% brightest pixels..

[counts, bins] = histcounts(sumIm(:), 256, 'Normalization', 'cumcount'); %this is CUMULATIVE COUNTS histogram. Therefore each bin column reflects total count less than its particular value
topPercent = 0.96*prod(size(imAdjRGB,1:2)); %identify the 96% percentile and above intensities. we are going to dampen those attention whores until they behave!
cnts = counts-topPercent;
cnts(cnts<0) = 0;
level = find(cnts, 1, 'first'); %these three lines identify for us the index with which we can manipulate exclusively the top percentage pixel int. 
idx = sub2ind(size(imAdjRGB,1:2), find(sumIm > bins(level)));

sumIm(idx) = bins(level); 
sumIm = rescale(sumIm); % this is an the all critical step. by dropping out the highest intensity pix's and then rescaling, homogenously gets brighters
sumIm(idx) = sImFilt(idx);

myMask = imextendedmax(sumIm, 0.7); %0.3 to 0.85

myMask = imclose(myMask, se1);
myMask = imopen(myMask, strel(ones(7,7)));

id = 'images:bwfilt:tie';
warning('off',id)
% get just one blob
myMaskOne = bwareafilt(logical(myMask),1);
maskymask = im2double(myMask);

areaBig = sum(maskymask(:));

%percentage of the whole mask you've made that is the 1 largest continuous blob
ratio = sum(myMaskOne(:))/areaBig;

if ratio < 0.96          %maybe two blobs or three blobs is better than 1. 10+ suggests a much bigger problem!  
    counter = 2;
    while ratio < 0.96 && counter < 11
        maskymask = bwareafilt(myMask,counter);
        ratio = sum(maskymask(:))/areaBig;
        counter = counter+1;
    end
    se3 = strel('Octagon',12);
else
  se3 = strel('Octagon',3);  
end
warning('on',id)

Ierode = imerode(maskymask,se3);
Ie_reconstruct = imreconstruct(Ierode, maskymask);
Ir_dilate = imdilate(Ie_reconstruct,se3);

Iopenbyrecon_closingbyrecon = imreconstruct(imcomplement(Ir_dilate),imcomplement(Ie_reconstruct));
Iopenbyrecon_closingbyrecon = imcomplement(Iopenbyrecon_closingbyrecon);

myMask = imregionalmax(Iopenbyrecon_closingbyrecon);

se2 = strel('diamond', 7);
binaryMask = imopen(myMask, se2);

if bweuler(binaryMask) < -10 %Euler number gets negative with more holes. So if its rather negative then we need to do an imfill
    binaryMask = imfill(binaryMask, 'holes');
end

end

