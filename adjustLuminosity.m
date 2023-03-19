function Afinal = adjustLuminosity(imgCurrent)

imgCurrent = im2single(imgCurrent);
matrxBlur = [0 -1 0; -1 5 -1; 0 -1 0;];
blurCurrent = imfilter(imgCurrent, matrxBlur, 'symmetric', 'conv');

%first make smaller because super slow
%img_small = imresize(imgCurrent, 0.75, {@oscResampling,4});

[L, alpha, beta] = fastRGB2Lab(blurCurrent, 0);
lum = L./100; lum = rescale(lum, 0, 1);
[B, ~, ~] = imreducehaze(lum,'ContrastEnhancement','boost', 'BoostAmount', 0.9);
[~, D] = imlocalbrighten(B,0.8, 'AlphaBlend',true);
Dinv = imcomplement(D);

LABim = cat(3, Dinv, alpha, beta);
[r,g,b] = fastLab2RGB(LABim);
adjRGB = (r+g+b)./3.0;


    sigma = 0.125; % bigger sigma reslts in duller, hazier images
    alpha = 6.0; %amount of smoothing increases with greater alpha. (do not set this to <1!!, that does opposite of smoothing)
    numLevels = 8; %fewer levels prioritizes speed over quality (range rec: [10, 100]
    beta = 6.0; % bigger geta means more dynamic range
    adjRGB = locallapfilt(adjRGB,sigma, alpha, beta,'NumIntensityLevels', numLevels);
    adjRGB = im2double(adjRGB);

adjRGB=imcomplement(adjRGB);
Afinal = imresize(adjRGB,size(imgCurrent,1:2));
        [level2, ~] = graythresh(ImageAdj);
        imBin = imbinarize(ImageAdj, level2);
        imBin2 = bwareafilt(imBin, 1);
        imB3 = imfill(imBin2, 'holes');
        outOfBoundsMask = activecontour(rgb2gray(imgCurrent), imB3, 15, 'edge', 'SmoothFactor' ,0.5);
    