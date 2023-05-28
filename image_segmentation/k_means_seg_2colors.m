function finalOutMask = k_means_seg_2colors(imgCurrent)
% outOfBoundsMask = k_means_seg_wrapperFcn(imgCurrent)
% This calls the image processing toolbox function IMSEGKMEANS,(a 
% segmentation method based on the k-means algorithm. It uses for this
% segmentation the RGB channels of the RGB colorspace, local entropy, 
% and local variance. we make use of the IMSEGKMEANS powerful abiliity 
% to handle any amount of channels, and that generally speaking,
% more information is only going to faciltate segmentation.
% We then use a (perhaps overkill) tristimulus vector distance method of identifying which
% color label is background and which is foreground. 

matrxBlur = [0 -1 0; -1 5 -1; 0 -1 0;]; %edge preserving (I think) kernel, blur fast.
rgbBlurred = imfilter(imgCurrent, matrxBlur, 'symmetric', 'conv');

Eim = rescale(im2single(cat(3, entropyfilt(imgCurrent(:,:,1)),entropyfilt(imgCurrent(:,:,2)),entropyfilt(imgCurrent(:,:,3)))),0,1);
Sim = im2single(cat(3, stdfilt(imgCurrent(:,:,1),ones(9)),stdfilt(imgCurrent(:,:,2),ones(9)),stdfilt(imgCurrent(:,:,3),ones(9))));
MedFilt = cat(3, medfilt2(imlocalbrighten(Eim(:,:,1)), [7, 7]), medfilt2(imlocalbrighten(Eim(:,:,2)), [7, 7]), medfilt2(imlocalbrighten(Eim(:,:,3)), [7, 7]));
    
chls = cat(3, rgbBlurred, Eim, rescale(Sim, 0,1), MedFilt);

numColors = 2;

LABELED = imsegkmeans(chls,numColors, 'Threshold', 1e-3); %the core of this function is the built in imsegkmeans
LABELEDrgb = im2double(label2rgb(LABELED));
outercolor = squeeze(mean(LABELEDrgb(2:8,2:5,:), [1,2])); 
    
nRows = size(imgCurrent,1); nCols = size(imgCurrent,2);
nPxls = prod(size(imgCurrent,1:2));
rgbStack = reshape(LABELEDrgb, nPxls, 3);

% Ensure that this is a row vector.
avgColor_vec = outercolor(:)'; 

% Compute the Euclidean distance between all rows of image and avgColor_vec. 
Dis = sqrt(sum((rgbStack - avgColor_vec(ones(nPxls,1),:)).^2, 2));

% Initialize I as a column vector.
Im = true(nPxls,1); 

% Set to 1 the locations in I at which D <= Thresh. 
Im(Dis <= 0.1) = false;

% Reshape I into an M-by-N image and fill in any holes
finalOutMask = imfill(reshape(Im,nRows,nCols), 'holes'); 

end
