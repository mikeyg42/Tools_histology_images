function outOfBoundsMask = runLevelSet(rgbImage)
% outOfBoundsMask = runLevelSet(rgbImage)
% This function implements a very simply formualtion of the level sets
% algorithm. syntax is: outOfBoundsMask = runLevelSet(rgbImage)
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% INPUT: an RGB Image, after preprocessing to correct for background ect.

% PROCESSING:
%  - image is significantly downsampled.
%  - Using gradients we extract a solid, boundary edge, which acts as a scalar
%    force field (when the expanding interface encounters white, it expands,
%    and when it encounters black it retreats).
%        - This is facilitated by some substantial pre-prcessing blurring).
%  - I initialize the zero-level set as a circle centered on the weighted centroid,
%    with a radius 1/6 the shorter side of the bounding box.
%  - The algorithm is allowed to cycle for 1250 iterations, which should be
%     overkill except instances with objects with long thin protrusions.
%  - For the sake of processing time, I opted against including a "reset"
%    as is recommended by the text of the reference text (see below). The reset would 
%    ensure the level-set function stays a sign-distance fuction throughout
%    iterations. But in testing my im's, segmentation worked without it well enough.
%  - After the level set algorithm completes, its obvious upon inspection that,
%    as a result of the blurring, we sometimes do initially to get a strong
%    boundary and as a result, our segmentation is not precisely lined-up with input
%    im. I correct this with a few iterations of a snake contour.
%  - Finally, i resize this mask back to the size of the input. 

% OUTPUT: a logical image of the same size as input, with
%       0 = background & 1 = foreground

%   NOTE: Sometimes this algorithm doesn't work... In that event, 
%   it will output a black rectangle with a white circle, which 
%   should prevent the whole Segementation GUI script to fail.
%
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% I prepared this implementation of Level Set Segmentation Algorithm by relying 
% on the following text:

% "Digital Image Processing Using Matlab" (4th ed., 2020) 
%   by Rafael Gonzalez, Richard Woods, and Steve Eddins. 
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

close all force

if ~isfloat(rgbImage)
    rgbImage = im2double(rgbImage);
end

BW = rgb2gray(rgbImage);

% this function is not a fast one... It helps a LOT though to work with a
% small image. I found it to run unbearably slow with array's of length > 1000
if max(size(BW, 1:2)) > 850
    BWsmall = imresize(BW, [750, NaN], {@oscResampling,4});
else
    BWsmall = BW;
end

BWsmall  = imadjust(BWsmall);
lg = fspecial('log', 7);
BWsmall = rescale(imadjust(imfilter(BWsmall, lg)), 0,1);

gaussianSmoothing = fspecial('gaussian', 15, 4);
fsmooth = imfilter(BWsmall, gaussianSmoothing, 'replicate');

% use the norm of the image gradient to identify the tissue's boundary edges.
% This boundary will become our scalar force field
[gy,gx] = gradient(fsmooth);
fnorm = sqrt(gx.^2 + gy.^2);

% arbitrarily set parameters fine-tuning the scalar field appearance
p = 1.5; lambda = 75; % we want lambda to be large enough to spread F's min and max to be maximally seperate

Fp = 1./(1 + lambda*(fnorm.^1.5)); %changed the 2 to 1.5????

%padarray so that later on, imfill works, but, also we want to pad outside of the loop, for speed
FprimePad = padarray(Fp, [3, 3], 1, 'both'); 

% define a threshold value for Fprime
T = (max(Fp(:))+min(Fp(:)))/2;

% this WHILE loop is an effort to try to make more robust/optimized the
% thresholding scheme for the scalar force field. Critically, the level set
% needs a boundary free of "leaking" or breaks. Thresh too low, lots of breaks. too high, its all black.
% I do a quick, albeit imperfect testing for gaps/leaks, namely to do a fill w/ IMFILL 
% -- if everything becomes white, there likely is a leak.
% -- if every value we check in the loop from T up to 1 has the same issue,
% then I've written a hail mary much more aggressive attempt at making F...

tf = 0;
k = 1.0;
while tf == 0
    F2 = FprimePad > T*k; % take the threshold 
    F2 = logical(F2);
    
    F_0 = imclose(imcomplement(F2), strel(ones(11))); %fill in any potential gaps in the boundary
    F_1 = bwmorph(F_0, 'thin', 2);
    F_1 = bwareaopen(F_1, 1500); %get rid of any objects with fewer than 1000 pixles
    
    forceField = imcomplement(F_1);
   
    test = imfill(forceField, 'holes'); 
    if sum(sum(test)) == prod(size(test, 1:2)) || sum(sum(test)) == 0 %if everything is zeros or ones, bad news
        k = k*1.02;
    else
        tf = 1; 
    end
    
    if T*k > 1 && tf ~=1 % we try one final time using the subfunction "tryAgainWithMoreContrast"
        forceField=tryAgainWithMoreContrast(BWsmall);   
        F_1 = imcomplement(forceField);

        test = imfill(logical(forceField), 'holes');
        if sum(sum(test)) == prod(size(test, 1:2)) || sum(sum(test)) == 0 %if it still fails... give up and return to segmentation menu
            disp('could not automatically threshold forcefield. returning fcn.');
            sz = round(size(BW, 1:2)./2);
            randomCircle = images.roi.Circle('Center',[sz(2), sz(1)],'Radius', min(sz)/3);
            outOfBoundsMask = createMask(randomCircle, BW);
            return
        else
            tf = 1;
        end
    end
end

% use the weighted centroid and the bounding box of the region to inform
% where you initialize the level set curve
     
% the try statement fails when size(BWsmall, 1:2) already equals size(F_1, 1:2), which
% occurs only if its necessary to run "tryAgainWithMoreContrast" where
% padding is lost
try props = regionprops('struct',F_1, padarray(BWsmall, [3,3], 1, 'both'), 'WeightedCentroid', 'BoundingBox');
    szFlag = 0;
catch 
    props = regionprops('struct',F_1, BWsmall, 'WeightedCentroid', 'BoundingBox'); 
    szFlag = 1; %this flag equalling one means size(forceField & F_1) is equal to BWsmall (and only BWsmalll
end
        
cPoint = props(1).WeightedCentroid;
x0 = round(cPoint(1)); y0 = round(cPoint(2)); %[x0, y0] is center of circle
bbox = props(1).BoundingBox;
BB_width = bbox(3); BB_height = bbox(4);
r = round(min(BB_width, BB_height)/6); % r is radius

% initallize phi by making a circle as described by [x0, y0] and r.
[m,n] = size(BWsmall, 1:2);
if szFlag == 0
    y = 1:n+6; x = (1:m+6)'; % accounting for the added padding!!!!
else %if 
    y=1:n; x = (1:m)';
end
phi0 = hypot(x - x0, y - y0) - r^2; 

% initiallize
phi = phi0;

% in lieu of passing the whole force field through this massive loop, we
% can input only what we need, which is these three values: 
minF = min(forceField, 0);
maxF = max(forceField, 0);
delT = 0.5*(1/max(abs(forceField(:))));
sz_padded_in_loop = size(phi, 1:2)+2;

% TIME FOR THE LOOP!
for I = 1:1200
    phi = levelsetIterate(phi, minF, maxF, delT, sz_padded_in_loop);
end

% this should rescale all values to range of [0, 1]
phi_scaled = adapthisteq(rescale(phi, 0, 1), 'clipLimit', 0.1,'Distribution', 'exponential');
phi_scaled = imcomplement(imtophat(imcomplement(phi_scaled), ones(9)));

% there are a lot of values in the e-17 range which just can be set to 0
phi_scaled(phi_scaled< 1e-5) = 0;

% When everything is working, the histogram of scaled phi always has a peak at 0,
% and then it also always has another peak slightly to the right of the zero peak,
% sometimes smaller, usually at a gray value around 0.2 or 0.3. We want to threshold right above that second peak!
phi_removingzeros = phi_scaled(:);
phi_nz = phi_removingzeros(phi_removingzeros>0);
Thresh = mode(phi_nz, 'all')*1.025;

if Thresh<0.45
edgeImage = phi_scaled > Thresh;
else
T = adaptthresh(adapthisteq(phi_scaled), 'ForegroundPolarity', 'dark','Statistic', 'gaussian', 'NeighborhoodSize', 2*floor(size(BWsmall)/32)+1);
edgeImage = imbinarize(adapthisteq(phi_scaled), T );
end
edgeImage = imclose(edgeImage, strel(ones(3)));
edgeImage = imopen(edgeImage, strel(ones(9)));

try finalMask = LinkUpEdgeDiscontinuities(imcomplement(edgeImage));
catch ME %if the mask has no discontinuities then it will fail and enter catch block
    if strcmp(ME.identifier, 'MATLAB:badsubscript') 
        finalMask = imcomplement(edgeImage);
    else
        rethrow(ME);
    end
end

% edgeImage looks like we've drawn a thick outline of the outside perimeter.
% We don't want that, we want a full mask, so we fill in the edge to make a solid mask.
mask = imfill(finalMask, 'holes');

%remove padding if necessary!
if szFlag == 0 
    mask2 =  mask(4:end-3, 4:end-3);
else
     mask2 =  mask;
end

% because we smoothed fairly aggressively, this mask won't necessarily
% coincide precisely with that of the unsmoothed image. So, we iterate a few times with
% a snake, refining the segmentation with now the unsmoothed image.
mask3 = activecontour(BWsmall, mask2, 36);

% note: if we resized before applying the active contour snakes instead of now, 
% the snakes processing time would go from less than a second to ~5min!
outOfBoundsMask = imresize(mask3, size(BW, 1:2));

end

%----------------------------------------------------------------------%

function phin = levelsetIterate(phi_nopad, minF, maxF, delT, sz)
%levelsetIterate caluclates the iterative solution of the level set equation.
%   Our our means of solving level set equation is iteratively. which we do here
% PHIN is the new iteration after PHI. Because we've chosen F (force term)
% invariant throughout iterations, and because we don't need all of F, only 3 
% descriptors from it, we save memory and input only PHI and those 3
% descriptors, which we've precalculated outside the loop. They are:
%       - delT, the time increment, given as 0.5*(1/max(F(:)) (should be
%        greater than 0 and less than or equal to 1.
%       - minF contains the negative values of F (positive values becoming 0), and
%           maxF contains only the positive values (the neg values becoming 0),
%           ie max(minF(:)) = min(max(F(:)) = 0.
%   syntax:  phin = levelsetIterate(phi, minF, maxF, delT)
%

%% COMPUTE UPWIND DERIVATIVES (EQ. (12-58) and THE UPWIND GRADIENT MAGNITUDES (EQ. (12-57)
%% step 1: upwind derivatives
% pad phi w/ 1 pixel so that you can calculate derivatives using the discrete numerical derivative
% method of central differences

phi = repmat(1, sz(1), sz(2));
phi(2:end-1, 2:end-1) = phi_nopad;

% vectorized way of calculating the 4 finite central differenes in a 4-way
% neighborhood around a point.
Dxplus  = circshift(phi, 1, 1) - phi;
Dxminus = phi - circshift(phi, -1, 1);
Dyplus  = circshift(phi, 1, 2) - phi;
Dyminus = phi - circshift(phi, -1, 2);

% Strip out the border
Dxplus  = Dxplus(2:end-1, 2:end-1);
Dxminus = Dxminus(2:end-1, 2:end-1);
Dyplus  = Dyplus(2:end-1, 2:end-1);
Dyminus = Dyminus(2:end-1, 2:end-1);
phi_rmPad = phi(2:end-1, 2:end-1);

%% step 2: upwind gradient magnitudes
%   There are two components of the upwind normalized gradient of a
%   level set function (+ and -). We can calcualte with the upwind
%   derivatives: Dxplus, Dxminus, Dyplus, and Dyminus from above.
gMagPlus  = sqrt((max(Dxminus,0).^2) + (min(Dxplus,0).^2) ...
    + (max(Dyminus,0).^2) + (min(Dyplus,0).^2));

gMagMinus = sqrt((max(Dxplus,0).^2) + (min(Dxminus,0).^2) ...
    + (max(Dyplus,0).^2) + (min(Dyminus,0).^2));

%% step 3: update Phi

delta = delT*(maxF.*gMagPlus + minF.*gMagMinus);

phin = phi_rmPad - delta;

end


function forceField=tryAgainWithMoreContrast(BWsmall)

bwEq = adapthisteq(BWsmall, 'Distribution', 'exponential'); % we revert to the unpadded!
fsmooth = imgaussfilt(imadjust(bwEq, [0.3, 0.7], [0, 1]), 1.25);

[gy,gx] = gradient(fsmooth);
Fprime = 1./(1 +  60*sqrt(gx.^2 + gy.^2)); %p = 1, lambda = 60

FprimePad = padarray(Fprime, [25, 25], 1, 'both'); 

T2 = (min(FprimePad(:))+max(FprimePad(:)))/2;
F2 = FprimePad > T2;
F3 = imclose(imcomplement(logical(F2)), strel(ones(11))); %fill in any potential gaps in the boundary

ff =  bwareaopen( bwmorph(F3, 'thin', 1), 1800);
FField = bwareafilt(imbinarize(imdilate(bwperim(imcomplement(ff), 8), ones(9)).*ff), 1);

forceField = FField(26:end-25, 26:end-25); %remove padding!

end
