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

imgPre = preprocessRawRGBims(imgIn, 0.5);
img = locallapfilt(im2single(imgPre), 0.02, 1.8,'NumIntensityLevels', 200);

se1 = strel(ones(3));
se2 = strel([0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0]);
se3 = strel('disk', 7);
% Nice way to increase contrast quickly
imgC = imsubtract(imadd(img,imtophat(img,se2)),imbothat(img,se2));
imgC = im2double(imgC);

%- BW1 ----------------------------------
tic;
lumIm = rgb2lab(imgC);
lumIm = double(lumIm(:,:,1)./100);

for k = 1:3
    Im = imcomplement(imgC(:,:, k));
    [gx, gy] = derivative7(Im, 'x', 'y');
    
    % Compute the magnitude image using hypot
    Gmag = imadjust(hypot(gx, gy));
    
    se11 = strel('disk',13);
    Ierode = imopen(imerode(Gmag,se11), se1);
    Ie_reconstruct = imreconstruct(Ierode, Gmag);
    Ir_dilate = imdilate(Ie_reconstruct,se11);
    Iopenbyrecon_closingbyrecon = imreconstruct(imcomplement(Ir_dilate),imcomplement(Ie_reconstruct));
    Gmag_recon = imclearborder(rescale(imcomplement(Iopenbyrecon_closingbyrecon), 0, 1));
    
    [n, edges] = histcounts(Gmag_recon, 0:0.05:1);
    listN = sortrows(vertcat(n, edges(2:end))',1, 'descend');
    bw1_thresh = (max(listN(1:2, 2))-0.05)-0.08;
    BW1a = Gmag_recon > bw1_thresh;
    
    %figure; imshowpair(Gmag_recon, BW1a, 'montage');
     temp = imfill(BW1a, 'holes');
     if (numel(find(temp(:)))-prod(size(BW1a, 1:2)))==0
          cBW1(:,:, k) = imclose(BW1a, strel('diamond', 7));
     else
          cBW1(:,:, k) = temp;
     end

end
cBW1 = im2double(cBW1);
BW1 = sum(cBW1,3)./3 == 1;

toc;
%%

%---- BW2 -------------------------------

tic; BW2 = fastMarchingEntropyFilter(lumIm); toc;

%------- BW3, 4, 5-----------------------
% texture filters: entropy, range, variance filters 
tic;
imsBinarized = threeTexturesSegmentation(imgC);

BW3 = imsBinarized{1};
BW4 = imsBinarized{2};
BW5 = imsBinarized{3};

toc;
%-------------- BW6 ---------------------
tic;
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
    channelGrayIms(:,:,k) = rescale(imSharp, 0, 1);
end
binaryUnion = channelBinaries(:,:,1) & channelBinaries(:,:,2) & channelBinaries(:,:,3);
BW6a = imfill(binaryUnion, 'holes')>0.8;
BW6b = hysteresisThreshold_wMorph(channelGrayIms);
BW6 = BW6a & BW6b; % this accomodates the fact that the BW6a might fail its imfill and fill the whole image with 1's.
figure; montage({BW6a,BW6b, BW6});
toc;
%----------------------------- BW9 ------
tic;
% this mess starts off using a classical background filter technique as others above did:
% (im+imtophat(im))-imbothat(im); this is followed by background correction by subtracting a huge
% median filter. 
% see: An efficient algorithm for retinal blood vessel segmentation using h-maxima transform and multilevel thresholding
%Saleh and Eswaran 2011 (this is where I got the idea, although the math itself is a century older)

imgA = lumIm+0.5.*imtophat(lumIm, se2);
imgD = imgA-0.5.*imbothat(imgA, se2);
imgB = imcomplement(histeqfloat(imgD - medfilt2(imgD, [25, 25])));

% Need to be careful about taking the log of an image because of log(0) issues creating
% infinities, also negatives inadvertently creating complex numbers. 
imgB = imclose(imgB, strel(ones(3)));
imgE = real(log( imgB));
imgE(imgE==Inf|imgE==-Inf) = eps; 
imgE = rescale(imgE, 0, 1);

% n=12 was determined through trial and error... may require some personal adjustment 
imgE = imbinarize(imgE.^12+hypot(derivative7(imgE, 'x'), derivative7(imgE, 'y')));
BW7 = imopen(imfill(imclose(imgE, se1), 'holes'), strel(ones(5))); 
BW7 = cleanMask(BW7);
toc;
%% PART 2: images are made = time to select the winners and combine their wisdom
%
ImList = {imgC, BW1, BW2, BW3, BW4, BW5, BW6, BW7};
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

pb1= uibutton(gl3_images,'Text', 'Confirm all good masks are checked, so we can move on',...
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
% rb8 = uicheckbox(gl2_buttons, ...
%     'Text', '     Opt8');
% rb9 = uicheckbox(gl2_buttons, ...
%     'Text', '     Opt9');

% layout of gl3 (the top panel)
ax.Layout.Row = [1 6];
pb1.Layout.Row = 7;
%layout of gl2 (bottom panel)
bg.Layout.Row = 5;
pnl.Layout.Row = [1, 4];

imshow(out, 'Parent', ax, 'Border', 'tight');
drawnow expose;

uiwait;

%query the checkboxes now!
hCheckboxes = findobj(gl2_buttons,'Type','uicheckbox');
checkboxValues = get(hCheckboxes, 'Value');

close all force

selectedImages = cell2mat(checkboxValues);
bestIms = ImList(selectedImages);
andImage = double(bestIms{1});
nIms = numel(bestIms);
for pp = 2:nIms
    andImage = andImage+double(bestIms{pp});
end

andImage(andImage<nIms)=0;
andImage(andImage==nIms)=1;
andImage = logical(andImage);
gg = regionprops('table', andImage, 'MajorAxisLength');
nRegions = sum(gg.MajorAxisLength>200);
finalMask = bwareafilt(andImage, nRegions);
warning('on','MATLAB:polyshape:repairedBySimplify');
end

function buttonCallback(~, ~)

uiresume;
end

function BW2 = fastMarchingEntropyFilter(lumIm)

%Use entropy filter. filter it to identify blobs of + signal. use pdist to find how far
%away each blob is from nearest neighbor, then apply outlier detection alg to that vector
%(I used >1.25 IQR). imreconstruct those blobs, then fit a tight boundary. Convert to
%polygon, and use this polygon as your seed for the fast marching algorithm.
lumIm = rescale(imboxfilt(imadjust(lumIm), 3, 'NormalizationFactor',1),0,1);

% Peter Kovesi's 7-stencil derivative
[gx2, gy2] = derivative7(rescale(entropyfilt(lumIm , true(7)), 0, 1), 'x', 'y'); 
Gmag = imlocalbrighten(hypot(gx2, gy2)); %hypot gives you the magnitude of your derivatives here 

Gfilt = Gmag-medfilt2(rescale(entropyfilt(lumIm , true(7)), 0, 1), [25, 25]);
Gfilt = imcomplement(imadjust(rescale(Gfilt, 0, 1)));

level = graythresh(Gfilt);
mask = imopen(imbinarize(Gfilt,level*1.25), strel('disk', 4,8));

GmagBlobs2 = imreconstruct(im2double(mask), imsharpen(Gfilt));
level = graythresh(GmagBlobs2);
maskBetter = imbinarize(GmagBlobs2,level*1.1);

%use a HUGE erosion to set your seed is definitely well within the blob, 
seed = bwareafilt(imerode(maskBetter, strel('disk', 80, 8)),1);

% anisotropic filtering is another edge preserving filter. histeqfloat is Peter Kovesi's
[fmm_im, ~] = imsegfmm(graydiffweight(histeqfloat(Gfilt), seed, 'GrayDifferenceCutoff', 25,  'RolloffFactor' , 0.4), seed, 0.005); %I just used the thresh and graydist cutoff that they do in documentation. 
fmm_processed = imfill(imreconstruct(maskBetter,fmm_im), 'holes'); 
fmm_processed = imreconstruct(im2double(fmm_processed), Gfilt);

centroid = regionprops('table',seed, 'centroid');
centroidSeed = round(centroid{1, :}); % recall that regionprops output is [X,Y]
graycon = grayconnected(fmm_processed, centroidSeed(2), centroidSeed(1), 0.25); % This functions wants rows/cols

BW2 = imresize(imfill(graycon, 'holes'), size(lumIm), 'bilinear');
end

function imsBinarized2 = threeTexturesSegmentation(imgC)

se2 = strel([0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0]);
se3 = strel('disk', 7);

Eim = imaverageChan(rescale(entropyfilt(imgC)));
Sim = imaverageChan(rescale(stdfilt(imgC,ones(13))));
Rim = imaverageChan(rescale(rangefilt(imgC,ones(13))));
imsTexture = {Eim, Sim, Rim};
imsBinarized = cellfun(@(x) imbinarize(x, 'adaptive'), imsTexture, 'UniformOutput', false);
imsBinarized = cellfun(@(x) imclearborder(x), imsBinarized, 'UniformOutput', false);
imsBinarized = cellfun(@(x) imclose(x, se2), imsBinarized, 'UniformOutput', false);
imsBinarized = cellfun(@(x) imfill(x,'holes'), imsBinarized, 'UniformOutput', false);
idx = false(1,3); nPixels = prod(size(imgC, 1:2));
for p = 1:3
    blobCC = bwconncomp(imsBinarized{p});
    cellArrayBlobs = cellfun(@length, blobCC.PixelIdxList');
        propMask = sum(cellArrayBlobs(:))/nPixels;
        if propMask<0.6 % 40% of the image at leeast should be real
            idx(1, p) = true; 
        end
        
        if size(cellArrayBlobs, 1)>=4
        bigblobs = maxk(cellArrayBlobs,4);
        idx(1, p) = mean(bigblobs(2:4)./bigblobs(1)) > 0.01;
        end
end
 % I take the size of the largest blob and divide the next 3 largest blobs with it, then average the proportions....
% If they are more than 1% of the big blob, that means im fill missed majority of it\
imsBinarized2 = imsBinarized; %preallocate with a replica
myFUNKS = {@(x) rescale(entropyfilt(x, ones(7)), 0, 1), @(x) rescale(stdfilt(x,ones(7)), 0,1), @(x)rescale(rangefilt(x,ones(7)),0,1)};
for k = 1:3
    if idx(1, k)
        imA = imsTexture{k};
        imB = imcomplement(im2single(imadjust(imA)));
        h = multithresh(imB, 7); 
        bw = imB<=h(1); imC = imcomplement(imimposemin(imB, bw));
        imC(imC>1e4) = median(imC(:));
        imD = locallapfilt(rescale(imC, 0, 1), 0.2, 3);
        myFunk = myFUNKS{k};
        BWfinal= cleanMask(imfill(imclose(imclearborder(imbinarize(myFunk(imD), 'adaptive')), se3), 'holes'));
        imsBinarized{k} = logical(BWfinal); 
    else
       imsBinarized2{k} = logical(cleanMask(imsBinarized{k}));
    end
end 
imsBinarized2 = cellfun(@(x) imfill(x, 'holes'), imsBinarized2, 'UniformOutput', false);
imsBinarized2 = cellfun(@(x) bwareafilt(x, 1), imsBinarized2, 'UniformOutput', false);
end

function img2 = imaverageChan(rgb)
[a,b,c] = imsplit(rgb);
img2 = rescale(imadd(imadd(a,b), c), 0, 1);
end