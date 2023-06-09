function  outputMask = vectorCalcSegment(inputImageRGB)
% outputMask = vectorCalcSegment(inputImageRGB)
%
% This is a segmentation function. It leerages 4 different discrete methods for numeric
% approximations of image derivatives or various orders. Using the contrast-enhanced luminosity
% image, we find approximations of the 0.6th, 1st, 2nd, and 2.6th order graidents. 
% We convert each of these into binary masks, and then we average them togerther before 
% thresholding the resulting composite image. After filling in this blob we have our final output 
% BW MASK!
%
% In a 5000x5000x3 image generating the mask takes 30-45 seconds. 
% 
% - Firstly, the most straightforward, this technique is using
% the built-in function imgradientxy to generate the first order partial derivatives.
% We can use these partial derivates (1 in regards to x and the other y) and calculate
% from them the overall 1st order gradient magnitude, using hypot (which handles overflow
% well).

% - Second, we use a likely more accurate estimation of the image gradient for 1st and 2nd
% order, using an implementation written by Peter Kovesi of the 7-tap approximation of the
% derivative. We combine the images produced by this approximation in order to generate
% the fsecond derivative in the gradient direction (SDGD) image, which is well known for
% strongely accentuating edeges. I take this oine step further and use image
% reconstruction with the local variance in order to emphasize not just any edges, but in
% particular edges with higher amounts of local variance, indicative of "real" edges. 

% - Thirdly, I have calculated two fractional derivatives using the numerical method of 
% Riemann and Liouville which was improved upon by the Caputo-Grümwald formula. Fractional
% Derivates often have advantageous edges and less halo-effects. They also offer us
% diversity by accentuarting different edges than integer derivaties. For my test images
% here I've found that the 0.6th and 2.6th orders achieve great results with my images,
% based on trial and error of values close by to x.5/2.5, which I reasoned would 
% offer the most novel info given what we aready have calculated.  

% Because derivative images are invariably noisey, care had to be taken to apply 
% edge-preserving blurring techniques alongside contrast enhancement. It may be that in 
% a different dataset these filter sizes and intensities might need to be shuffled around.  

% NOTE: I tried to as much as possible avoid using any morphological dilation/erosion, so
% to keep the mask as sharp and precise as possible. I have permitted just one, for the
% SDGD processing, as nothing else made sense. 

% Michael Glendinning, May 2023
%
% ======================================================================== %

%% preprocess
labim = rgb2lab(inputImageRGB);
grayIm = histeqfloat(rescale(labim(:,:,1), 0, 1));

%% SDGD
% measure local variance for SDGD enhancement
localVarianceMap = medfilt2(stdfilt(grayIm, true(7)));
lVarMap = imadjust(imhmin(localVarianceMap, 0.15)); %supresses relative min less than 0.15 - reduces noise! 
T = multithresh(lVarMap, 5);

% Calculate the SDGD - Second Derivative in the gradient direction
[dx, dy, dxy, dxx, dyy] = derivative7(grayIm, 'x', 'y', 'xy','xx', 'yy');
SGDG = ordfilt2((dxx.*dx.^2+dyy.*dy.^2+2.*dxy.*dx.*dy)./(dx.^2+dy.^2), 9, ones(3)); %maximum filter to counteract DULLNESS

% use morphological reconstruction to prioritize edges in areas with high variance
edgeMap = imreconstruct(edge(imclose(rescale(imadjust(SGDG), 0,1), strel(ones(3))), 'Sobel'), imbinarize(lVarMap , T(3)));
im3 = bwmorph(bwareafilt(edgeMap, 1), 'thin');
if checkImfill(im3)
    im3 = imfill(im3, 'holes');
    if bwarea(im3)<0.5*prod(size(im3, 1:2))
        [y, x] = find(bwskel(im3));
        cent = round(sum([x, y], 1)./length(x));
        im3 = grayconnected( im2double(im3),cent(1), cent(2), 0.5);
    end
end

%% LOW-order fractional deriv (v = 0.6)
fractDeriv1 = imFractionDeriv(grayIm, 0.6);
fDeriv1 = imcomplement(histeqfloat(locallapfilt(im2single(fractDeriv1), 0.3, 0.3)));
T1 = graythresh(fDeriv1)*1.2;
im1 = imbinarize(fDeriv1, T1);

%% HIGH-order fractional deriv (v = 2.6)
fractDeriv2 = imFractionDeriv(grayIm, 2.6);
fDeriv2 = imcomplement(histeqfloat(medfilt2(fractDeriv2))); 
T2 = graythresh(fDeriv2)*1.2;
im2 = imbinarize(fDeriv2, T2);

%% Classic, built-in IPT gradient image (1st order)
[gx, gy] = imgradientxy(imlocalbrighten(grayIm));
gmag = medfilt2(hypot(gx, gy), [5,5]);
im4 = imbinarize(gmag, graythresh(gmag)*0.95);
im4 = bwareafilt(im4,1);
if checkImfill(im4)
    im4 = imfill(im4, 'holes');
end

%% Composite image formation
sumIM = im2double(im1)+im2double(im2)+im2double(im3)+im2double(im4);
sumIM_bin = imbinarize(rescale(sumIM,0,1), 0.4);
outputMask = bwareafilt(imfill(sumIM_bin,'holes'), 1);

end