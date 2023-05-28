function [BWmask, BWedge] = fullColorSegmentation(rawRGBim)
% BWmask = colorVectorSegment(RGBim)
% This function segments RGB images with white background. My goal here was to capture
% edges in color that might be missed by calling RGB2GRAY or any similar colorspace
% intensity projection. To start, we begin by seperating the RGB channels, increasing the
% contrasts of each, padding them, and then applying a sliding correlational spatial
% filter: namely, we collect a variant of the typical high pass filter by taking a sliding
% 5x5 neighborhood and calculatig for each pixel the deltaE magnitude of difference of
% each tristumulus from the average of the 24 pixels around it. We then take the
% derivative of this deltaE image, using Peter Kovesi's funciton for the 7-tap estimate of
% image partial derivaties. We apply local contrasting via IMLOCALBRIGHTEN to the magnitude of the
% second derivattive gradient.

% After this, we use image reconstruction to enhance our 1st derivative image with the
% enhanced edges of the second derivative. After this, using the morphological function
% IMIMNPOSEMIN applied to this enhanced image achieves an image very amenable to
% segmentation!

% Michael Glendinning, Apr 2023
% ============================================================================================

% Extract each color channel
[r, g, b] = imsplit(im2single(rawRGBim));
channels = {r, g, b};

% Compute the highpassfiltering using in a 5x5 neighborhood using colfilt
adjusted_channels = cellfun(@imadjust, channels, 'UniformOutput',false);
pad_channels = cellfun(@(x) padarray(locallapfilt(x, 0.8, 1.25), [9,9], 'replicate', 'both'), adjusted_channels, 'UniformOutput',false);
neighborhood_sums = cellfun(@(x) colfilt(x, [5, 5], 'sliding', @sum), pad_channels, 'UniformOutput',false);
highpassfilt = arrayfun(@(x) (neighborhood_sums{x} - pad_channels{x})./24, 1:3, 'UniformOutput',false);

% Use deltaE approximation to quantify relative differences in color
dE = rescale(sqrt(highpassfilt{1}.^2 + highpassfilt{2}.^2 + highpassfilt{3}.^2), 0,1);

% Take derivative to locate edge magnitude, then remove padding
[ex, ey] = derivative7(dE, 'x', 'y');
edge_mag = hypot(ex, ey);
contrast_edgeMag = imlocalbrighten(rescale(edge_mag, 0, 1));

%morphological reconstruction of the 1st deriv with the 2nd deriv
highCon = imreconstruct(contrast_edgeMag, dE) ;

% Using the mean of the image as basis for unsupervised im processing 
id_mean = find(round(highCon, 4) == round(mean(highCon(:)),4));
meanPoints = zeros(size(highCon, 1:2), 'single');
meanPoints(id_mean) = 1;

% morphological imimposemin! then unsupervised thresholding using the image mode which
% will be the background usually! 
flatim = imimposemin(imcomplement(highCon), meanPoints);
noPadding = flatim(10:end-9, 10:end-9); %remove padding
BWmaskp = noPadding > mode(noPadding(:))*1.01;
BWmaskp = bwareaopen(BWmaskp, 120);

BWmask_sing = activecontour(rawRGBim, BWmaskp, 12, 'Chan-Vese');

BWmask= im2double(BWmask_sing);

BWedge = bwperim(BWmask);

end