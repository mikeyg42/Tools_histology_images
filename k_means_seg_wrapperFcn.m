function outMask = k_means_seg_wrapperFcn(imgCurrent,varargin)
% outOfBoundsMask = k_means_seg_wrapperFcn(imgCurrent, [y/n to incl textures])
% This calls the image processing toolbox function IMSEGKMEANS,(a 
% segmentation method based on the k-means algorithm. It uses for this
% segmentation the *a and *b channels of the L*a*b colorspace (which
% contain the images color information, sans luminence.

% If 1 or true is included as an input to this function, then some
% additional channels are appended to this image stack, making use of the
% IMSEGKMEANS useful abiliity to handle any amount of channels, given that
% more information is only going to faciltate segmentaion. It uses STDFILT,
% which is the local standard deviation in a 9x9 nhood, and ENTROPYfilt,
% which caluclates the local entropy within that same size nhood. It helps
% considerably and in most cases should be included


try tf = islogical(varargin{1});
   if tf
    texture_y_n = varargin{1}; end
catch
    texture_y_n = 0;
end

matrxBlur = [0 -1 0; -1 5 -1; 0 -1 0;]; %edge preserving (I think) kernel, blur fast.
rgbBlurred = imfilter(imgCurrent, matrxBlur, 'symmetric', 'conv');

labCurrent= rgb2lab(rgbBlurred);
ab1  = labCurrent(:,:,2)  ; ab2  = labCurrent(:,:,3); %ignore luminosity, just focusing on color
chl1 = im2single(ab1)     ; chl2 = im2single(ab2)   ;

if ceil(mean(chl1(2:5,2:5), 'all'))==1 %check in the top left corner and see if its more or less white or black. flip if white
    chl1 = imcomplement(chl1); end
if ceil(mean(chl2(2:5,2:5), 'all'))==1
    chl2 = imcomplement(chl2); end

if texture_y_n == 1
    Eim = im2single(entropyfilt(imgCurrent));
    Sim = im2single(stdfilt(imgCurrent,ones(9)));
    chls = cat(3, chl1, chl2, Eim, Sim);
elseif texture_y_n == 0
    chls = cat(3, chl1, chl2);
end

numColors = 2;

LABELED = imsegkmeans(chls,numColors, 'Threshold', 1e-3); %the core of this function is the built in imsegkmeans

labelVals = [1, 2];
zerocolor = mean(LABELED(2:5,2:5), 'all');
maskColor = labelVals(labelVals~=zerocolor);

LABELED(LABELED==maskColor) = 5;
LABELED(LABELED==zerocolor) = 0;
maskBW = double(LABELED)./5;

se = strel('disk', 7, 8);

BW      = imopen(maskBW, se);
outMask = imfill(imcomplement(BW), 'holes');
outMask = logical(outMask);


