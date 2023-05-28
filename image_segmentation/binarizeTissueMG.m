function BinMask = binarizeTissueMG(RGBim, varargin)
%BinMask = binarizeTissueMG(RGBim) OR BinMask = binarizeTissueMG(RGBim, bufferAmount)
% this function will quickly binarize an RGB image. This is accomplished by extracting the
% luminosity from the image, flattening and blurring it with edge-preserving functions,
% followed by a generous thresholding. After this, we select just the largest blob,
% calculate its contour, which we turn into a polyshape. The polyshape allows us to take
% advantage of the polybuffer function, which I like more than imdilate for this. After
% we've applied the polybuffer, I convert it back to a mask, to which I apply 32
% iterations of an active contour snake, using as my base image the inverse of the
% gradient of the luminosity channel. 
% What this all amounts to is a relatively robust, relatively fast, uncomplicated means of
% segmenting a simple RGB image. It works best when the background is clear and flat. 

% Michael Glendinning, 2021

%% parsing optional input of polybuffer amount (default is 9)
try bufferAmount = varargin{1};
catch
    bufferAmount = 9;
end

%% relies on imfill, imdiffusefilt, and the roicolor on the luminosity chl (L*A*B image)

lab = rgb2lab(RGBim);
lscaled = rescale(lab(:,:,1), 0, 1);
lCh = imcomplement(lscaled);

% base the size of the kernel on the rough size of the image
if max(size(lCh, 1:2))> 6000
    se = strel('disk', 21, 8);
else
    se = strel('disk', 17, 8);
end

imGray = imfill(imclose(lCh, se));
% this imclose and fill helps bring areas with less luminance (such as a white matter lesion) up to the threshold for binarization,

% binarization of gray image using roicolor
[Tlevel, ~] = graythresh(imGray);
Idiffusion = imdiffusefilt(imGray);
newROImask = roicolor(Idiffusion, (Tlevel/1.5), 1);

% select the largest blob
BW = bwareafilt(logical(newROImask), 1);

%define the blob boundaries using contour
ct = contourc(im2double(BW), [0.5, 0.5]);
conX = ct(1,2:end);
conY = ct(2,2:end);

%convert contour into a polyshape and simplify the edges
ps = polyshape(conX, conY, 'Simplify', true);

%buffer the edges a bit to smoothen some jagged edges caused by the fat SE.
% then, remake your mask
ps2 = polybuffer(ps, bufferAmount);
bufferedMask = poly2mask(ps2.Vertices(:, 1), ps2.Vertices(:, 2), size(RGBim, 1), size(RGBim, 2));

%remove your padding of edges with 25 iterations of active contour
[gx, gy] = derivative7(lscaled, 'x', 'y');
target = imcomplement(hypot(gx, gy));
BinMask = activecontour(target, bufferedMask, 21);
end

