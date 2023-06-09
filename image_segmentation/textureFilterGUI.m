function finalMask = textureFilterGUI(imgIn, varargin)
% theMask = textureFilterGUI(img)
% A reliable semi-automated segmentation tool for grayscale images
% Grayscale images are input, and then 8 different possible binary
% masks are generated. Most of these segmentation techniques are classic
% in the field and implemented without substantial conceptual modifcation. In various ways
% I set parameters optimally for my data, so if poor performance in encountered, you will
% want to run through the code for those instances.
% There are 8 images to pick from:
%   - 1. image gradient, and is implemented based on Matlab documentation
%   - 2. some manipulation of the first image to
%   ensure that I've segemented the largest 2.5% of blobs in the image.
%   - 3, 4, 5. all rely on different means of locally filtering
% in the spatial domain so to assess regional differences in specific textures.
% I've kept the convolutions here at the recommended 9x9 square kernel. I
% have entropy(areas exhibiting randomness are brighter), range (brighter =
% more range surrounding a particular pixel), and std of course has
% brighter pixels indicative of more variation in neighborhood
% .  - 6.  I have the very classical gray-level co-occurance matrix
% derived image
%
%    - 7. 8. imextendedmax and min use the distance transform between known
%    boundary points and then converts to H-maxima transforms.
%    similar to the watershed seg. algorigthm

% Michael Glendinning, 2023
warning('off','MATLAB:polyshape:repairedBySimplify');

% filter with edge aware bilateral filt.
img = locallapfilt(im2single(imgIn), 0.02, 1.8,'NumIntensityLevels', 225);

% predefine these for later
se1 = strel(ones(3));
se2 = strel([0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0]);
se3 = strel('disk', 7);
se11 = strel('disk',11);

% Nice way to increase contrast quickly
imgC = imsubtract(imadd(img,imtophat(img,se2)),imbothat(img,se2));
imgC = padarray(im2double(imgC), [7,7],'replicate', 'both');
szPad = size(imgC, 1:3);
%- BW1 ----------------------------------
ImR = imgC(:,:, 1);ImG = imgC(:,:, 2);ImB = imgC(:,:, 3);

Im = cat(2, ImR, ImG, ImB);
    %calculte the first partial derivative with respect to x and y seperately
    [gx, gy] = derivative7(Im, 'x', 'y');
    
    % Compute the 1st deriv magnitude image using hypot and gx, gy
    Gmag = imadjust(hypot(gx, gy));
    
    Ierode = imopen(imerode(Gmag,se3), se1);
    Ie_reconstruct = imreconstruct(Ierode, Gmag);
    Ir_dilate = imdilate(Ie_reconstruct,se11);
    chlsMontage = imreconstruct(imcomplement(Ir_dilate),imcomplement(Ie_reconstruct));
    
ncol = szPad(2);
RGBadj_pad = zeros(szPad, 'double');
RGBadj_pad(:,:,1) = rescale(chlsMontage(:, 1:ncol), 0,1);
RGBadj_pad(:,:,2) = rescale(chlsMontage(:, 1+ncol:2*ncol), 0,1);
RGBadj_pad(:,:,3) = rescale(chlsMontage(:, 1+2*ncol:end), 0,1);

[BW1b, binEdge] = fullColorSegmentation(RGBadj_pad);

% check for if you need to invert image
test = imclearborder(BW1b);
if test(1:20, 1:20) ~= BW1b(1:20, 1:20)
BW1b = imcomplement(BW1b);
end

% check for leaky border before edge filling
test2 = imfill(BW1b, 'holes');
if test2(1:15, 1:15) == BW1b(1:15, 1:15)
    BW1 = test2;
else
    BW1 = BW1b | binEdge;
end
%%

%---- BW2 -------------------------------
chls = {ImR, ImG, ImB};
bwIms = cellfun(@(x) fastMarchingEntropyFilter(x),chls, 'UniformOutput', false);
[BW2a, successYN] = smartCombineChannels(bwIms);

if successYN 
    BW2 = BW2a;
else
    BW2a = fastMarchingEntropyFilter(rgb2gray(imgIn));
    BW2 = padarray(BW2a, [7, 7], 'replicate', 'both');
end

%------- BW3, 4, 5-----------------------
% texture filters: 3 = entropy, 4 = range, 5 = variance filters 
imsBinarized = threeTexturesSegmentation(imgC);

BW3 = imsBinarized{1};
BW4 = imsBinarized{2};
BW5 = imsBinarized{3};

%-------------- BW6 ---------------------
% We go through each channel of the RGB image. We get the  gray-level 
% co-occurance matrix on the average of the 1st and 2nd xy derivatives. 
% Then we blur slightly before a h-max transform. From this extract 2 images:
% 1. a grayconnected binary image (based on the centerpoint of the whole image)
%    After the loop we find the union of all three channel binarys
% 2. the grayscale image as is. After the loop we treat the image as an RGB im,
%    (which it is, albeit with distorted colors) and segment it with my 
%     hysteresis function, which takes advantage of similar processes as this fcn.

channelGrayIms = zeros(size(imgC));
channelBinaries = channelGrayIms;

sharpFilt = [0 0 0;0 1 0;0 0 0]-fspecial('laplacian',0.2);
% we will need this to specify the seed point for GRAYCONNECTED later.... 
szc = floor(size(imgC, 1:2)./2);
for k= 1:3
    [gx, gxx, gy, gyy] = derivative7(imgC(:,:,k), 'x', 'xx', 'y', 'yy');
    deriv1 = imadjust(rescale(hypot(gx, gy), 0,1));
    deriv2 = imadjust(rescale(hypot(gxx, gyy), 0,1));
    derivAvg = (deriv1+deriv2)./2;
    offsets = [0 1; -1 1;-1 0;-1 -1];
    [~,scaledI] = graycomatrix(derivAvg, 'Offset',offsets, 'Symmetric', true);
    gray = rescale(scaledI, 0, 1);
    im = imdiffusefilt(imclose(gray, se2), 'NumberOfIterations' , 20, 'ConductionMethod','quadratic');
    im = imgaussfilt(im, 1.2);
    im_max = imhmax(im, 0.9);
    im_max = im_max./max(im_max(:));
    imSharp = imfilter(im_max,sharpFilt,'symmetric');
    channelBinaries(:,:,k) = grayconnected(imSharp, szc(2), szc(1), 0.01);
end
BW6 = channelBinaries(:,:,1) | channelBinaries(:,:,2) | channelBinaries(:,:,3);

%----------------------------- BW7 ------
% this mess starts off using a classical background filter technique as others above did:
% (im+imtophat(im))-imbothat(im); this is followed by background correction by subtracting a huge
% median filter. 
% see: "An efficient algorithm for retinal blood vessel segmentation using h-maxima
% transform and multilevel thresholding" {Saleh and Eswaran 2011} 
% This article is what inspired me to proceed as I have, although ultimately only a small
% amount is the same. the strategy of course subserving the method is more than a century older.
lum = rgb2lab(padarray(imgIn, [7,7], 'replicate', 'both'));
lumIm = rescale(lum(:,:,1), 0, 1);
imgA = lumIm+0.5.*imtophat(lumIm, se2);
imgD = imgA-0.5.*imbothat(imgA, se2);
imgB = imcomplement(histeqfloat(imgD - medfilt2(imgD, [25, 25])));

% taking the log of an image demands some care so to avoid creating Inf's by taking log(0)
% and perhaps worse, creating complex numbers by trying to take log(negative numbers)!
imgB = imclose(imgB, strel(ones(3)));
imgE = real(log(imgB));
imgE(imgE==Inf|imgE==-Inf) = eps; 
imgE = rescale(imgE, 0, 1);

% n=12 was determined through trial and error... may require some personal adjustment 
imgE = imbinarize(imgE.^12+hypot(derivative7(imgE, 'x'), derivative7(imgE, 'y')));
BW7 = imopen(imfill(imclose(imgE, se1), 'holes'), strel(ones(5))); 
BW7 = cleanMask(BW7);

%% PART 2: images are made = time to select the winners and combine their wisdom
%
ImList = {BW1, BW2, BW3, BW4, BW5, BW6, BW7};
ImListpad = cellfun(@(x) padarray(x, [25, 25], 1),ImList, 'UniformOutput', false);
out = imtile(ImListpad, 'GridSize', [3, 3], 'BorderSize', 25);

figy = uifigure('Position', [350   150   600   700]);
gl_master = uigridlayout(figy, [5, 1]);
bg = uibuttongroup(gl_master);

pnl = uipanel(gl_master);
gl3_images = uigridlayout(pnl, [7, 1]);

ax = uiaxes(gl3_images);
ax.GridLineStyle = 'none';
disableDefaultInteractivity(ax);

pb1= uibutton(gl3_images,'Text', 'Once the best few options are checked, click here to combine them and proceede',...
    'ButtonPushedFcn',@buttonCallback);

gl2_buttons = uigridlayout(bg, [3, 3]);
rb1 = uicheckbox(gl2_buttons, ...
    'Text', '     Opt1');
rb2 = uicheckbox(gl2_buttons, ...
    'Text', '     Opt2');
rb3 = uicheckbox(gl2_buttons, ...
    'Text', '     Opt3');
rb4 = uicheckbox(gl2_buttons, ...
    'Text', '     Opt4');
rb5 = uicheckbox(gl2_buttons, ...
    'Text', '     Opt5');
rb6 = uicheckbox(gl2_buttons, ...
    'Text', '     Opt6');
rb7 = uicheckbox(gl2_buttons, ...
    'Text', '     Opt7');

% layout of gl3 (the top panel)
ax.Layout.Row = [1 6];
pb1.Layout.Row = 7;

%layout of gl2 (bottom panel)
bg.Layout.Row = 5;
pnl.Layout.Row = [1, 4];

imshow(out, 'Parent', ax, 'Border', 'tight');
drawnow expose;

uiwait;

%query the checkboxes !
hCheckboxes = findobj(gl2_buttons,'Type','uicheckbox');
checkboxValues = get(hCheckboxes, 'Value');

close all force

selectedImages = cell2mat(checkboxValues);
bestIms = ImList(selectedImages);
andImage = double(bestIms{1}); %initialize with the first bestIm (specifying double so they add)
nIms = numel(bestIms);
for pp = 2:nIms
    andImage = andImage+double(bestIms{pp});
end

andImage(andImage<nIms)=0;
andImage(andImage==nIms)=1;
andImage = logical(andImage);

cc = bwconncomp(andImage, 8);
sz = cellfun('prodofsize', cc.PixelIdxList)>400;
nRegions = numel(find(sz));
finalMask = bwareafilt(andImage, nRegions);
finalMask = finalMask(8:end-7, 8:end-7);

warning('on','MATLAB:polyshape:repairedBySimplify');
end

function buttonCallback(~, ~)
uiresume;
end

function BW = fastMarchingEntropyFilter(myImage)

% extract the luminosity and get local area sums via a mean filter (also called box
% filter) with normalization factor set to 1. 
myImage = rescale(imboxfilt(imadjust(myImage), 3, 'NormalizationFactor',1),0,1); %note that if lumIm is not floating, this will end up w bad overflow!!

% Find the edges w the 1st derivative, via Peter Kovesi's 7-stencil 1st part deriv code
[gx2, gy2] = derivative7(rescale(entropyfilt(myImage , true(7)), 0, 1), 'x', 'y'); 
Gmag = imlocalbrighten(hypot(gx2, gy2)); %hypot gives you the magnitude of your derivatives here 

% Calculate the entropy filter in a 7x7 neighborhood. Use a large median filter to edge-preserve blur
% and then combine this filtered entropy-image with the 1st derivative magnitude 
Gfilt = Gmag-medfilt2(rescale(entropyfilt(myImage , true(7)), 0, 1), [25, 25]);
Gfilt = imcomplement(imadjust(rescale(Gfilt, 0, 1)));

% Now in order to create a seed for the fast marching, we use morphological manipulations
% to first consolidate the blob (isolating the best one) and then eroding a LOT.
level = graythresh(Gfilt);
mask = imopen(imbinarize(Gfilt,level*1.25), strel('disk', 4,8));
GmagBlobs2 = imreconstruct(im2double(mask), imsharpen(Gfilt));
level = graythresh(GmagBlobs2);
maskBetter = imbinarize(GmagBlobs2,level*1.1);
seed = bwareafilt(imerode(maskBetter, strel('disk', 80, 8)),1);% HUGE erosion step ensure seeds is well within target blob, 

% Using the Gfilt composite image as a mask and the eroded blob as seed, perform fast
% marching algorithm!!
[fmm_im, ~] = imsegfmm(graydiffweight(histeqfloat(Gfilt), seed, 'GrayDifferenceCutoff', 25,  'RolloffFactor' , 0.4), seed, 0.005); %I just used the thresh and graydist cutoff that they do in documentation. 

% processing the algorithm's results
fmm_processed = imfill(imreconstruct(maskBetter,fmm_im), 'holes'); 
fmm_processed = imreconstruct(im2double(fmm_processed), Gfilt);

% use gray connected region texture aware filter to further enhance this result
centroid = regionprops('table',seed, 'centroid');
centroidSeed = round(centroid{1, :}); % recall that regionprops output is [X,Y]
graycon = grayconnected(fmm_processed, centroidSeed(2), centroidSeed(1), 0.25); % This functions wants rows/cols

BW = imresize(imfill(graycon, 'holes'), size(myImage), 'bilinear');
end

function imsBinarized2 = threeTexturesSegmentation(imgC)

%se2 = strel([0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0]);
se3 = strel('disk', 6);

Eim = imaverageChan(rescale(entropyfilt(imgC)));
Sim = imaverageChan(rescale(stdfilt(imgC,se3.Neighborhood)));
Rim = imaverageChan(rescale(rangefilt(imgC,se3.Neighborhood)));

[n, edges] = histcounts(Eim);
n = movmedian(n,3);
nInd = max(islocalmin(n));
thresh = (edges(nInd)+edges(nInd+1))/2;
bIm1 = imbinarize(Eim, thresh);

imsTexture = {Sim, Rim};
imB = cellfun(@(x) imcomplement(im2single(imadjust(x))), imsTexture, 'UniformOutput', false);
h1 = multithresh(imB{1}, 5); h2 = multithresh(imB{2}, 5); 
bw1 = imB{1}<=h1(2); bw2 = imB{2}<=h2(2); 

imgC1 = imcomplement(imimposemin(imB{1}, bw1));
imgC1(imgC1>1e4) = median(imgC1(:)); 
imD1 = locallapfilt(rescale(imgC1, 0, 1), 0.2, 3);

imgC2 = imcomplement(imimposemin(imB{2}, bw2));
imgC2(imgC2>1e4) = median(imgC2(:)); 
imD2 = locallapfilt(rescale(imgC2, 0, 1), 0.2, 3);

bIm2 = imbinarize(imD1);
bIm3 = imbinarize(imD2);
imsBinarized2 = {bIm1, bIm2, bIm3};

end


function img2 = imaverageChan(rgb)
[a,b,c] = imsplit(rgb);
img2 = rescale(imadd(imadd(a,b), c), 0, 1);
end


function [BWimage, flag] = smartCombineChannels(cellArray3Ims) 

% define a range of percentages between which each image fills the page
area100 = prod(size(cellArray3Ims{1}, 1:2));
rangeA = [area100*0.93, area100*0.67];
areaVals = cellfun(@(x) bwarea(x), cellArray3Ims, 'UniformOutput', true);

% index only those values within range
idx = areaVals<rangeA(1) & areaVals>rangeA(2);
successsYN = 1;

switch sum(idx)
    case 0
        BWimage = zeros(size(cellArray3Ims{1}));
        successsYN = 0;
        
    case 1    
        BWimage = cellArray3Ims{idx};
        
    case 2 
        ix = find(idx);
        BWimage = immultiply(cellArray3Ims{ix(1)}, cellArray3Ims{ix(2)}); %elementwise  multiplication

    case 3
        BWimage = cellArray3Ims{1} & cellArray3Ims{2} & cellArray3Ims{3};

    otherwise
        error('the ostensibly simple smartCombineChannels function is totally whack? this should not be possible!!!')
end

flag = logical(successsYN); 

end