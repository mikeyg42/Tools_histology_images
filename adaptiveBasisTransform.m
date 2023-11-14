function newImg = adaptiveBasisTransform()

ids = imageDatastore("/Users/mikeglendinning/projects/imageProcessing_github/test_images/");
imgCurrent = im2single(imresize(readimage(ids,3),0.8));

%background correction
imgCurrent = evenFasterCorrection(imgCurrent);

% the reason why the usual kmeans implementations have soe blurring is to create broader
% and more homogeneous swaths of color  - helps segment but not for this
matrxBlur = [0 -1 0; -1 5 -1; 0 -1 0]; %edge preserving (I think) kernel, blur fast.
ig = imfilter(imgCurrent, matrxBlur, 'symmetric', 'conv');

newTiny = imresize(imgCurrent, 0.25, "lanczos3");
newTiny = imresize(newTiny, [16, 16], "lanczos3");


%% focus on the background color first
newTiny(2:15, 2:15, :) = NaN;

% use entropy to identify those squares (out of 16x16 grid) that do not partially
% intersect the outer rim of grid     
Eim = mean(cat(3, entropyfilt(imgCurrent(:,:,1)),...
    entropyfilt(imgCurrent(:,:,2)),...
    entropyfilt(imgCurrent(:,:,3))), 3);
% [note, output of entropy filt is always double, even though imgCurrent is single]

% Resize same grid as RGB image,
EimEdge = imresize(Eim, [16, 16], "lanczos3");

% Set to nan the inner blocks of the entropy grid
EimEdge(2:15, 2:15, :) = NaN;

% Sort the entropy levels into quartiles, take the half lowest
eQuants = quantile(rmmissing(EimEdge(:)),0.4); % the lower the number, the more things that will be set to NaN
EimEdge(EimEdge>eQuants(1)) = NaN;

%% use the average of all the lowest entropy outer grid cells to get background estiation

% create a logical array to index into newTiny to extract those blocks withe the lowest
lowEntEdge_idx = ~isnan(EimEdge); %this logical array corresponds to the locations of low entropy
lowEntEdge_idx = repmat(lowEntEdge_idx, [1, 1, 3]);

% Use this index on the low entropy blocks of newTiny
edgeVals = newTiny(lowEntEdge_idx);
edgeValList = reshape(edgeVals, [],  3);

% taking the mean of this should give us a decent estimate of the background
edgeVal = mean(edgeValList, 1);

%% now we need to identify the other two colors (which we acutually care about)

sz = size(imgCurrent);
pxls = reshape(imgCurrent, prod(sz(1:2)), sz(3));

% This results in an OK tri-modal seperation
dis = pdist2(edgeVal, pxls, 'minkowski', 1.5);

% we want to separate out those values closest to the edge color.
thresh = findValley(dis);

% Remove all x-values that are = thresh or lower than thresh.
closeToEdgeIdx = dis <= thresh; % Pixels close to the edge color

% Get the remaining pixels that are not close to the edge color
remainingPxls = pxls(~closeToEdgeIdx, :);

% Now, remove the corresponding distances
dis_left = dis(~closeToEdgeIdx);

% Calculate histogram data
[counts, edges] = histcounts(dis_left, 'BinMethod', 'fd');
binCenters = edges(1:end-1) + diff(edges) / 2; % get the bin center

% Smooth the histogram counts data so that there are exactly 2 peaks
smoothedData = smoothSignal_MG(counts);
[~, locs] = findpeaks(smoothedData);

% Find the distances closest to color1idx and color2idx
[~, minIdx1] = min(abs(binCenters - binCenters(locs(1))));
[~, minIdx2] = min(abs(binCenters - binCenters(locs(2))));
%ddzzues closest to the bin centers
%fidentified by the peaks
closestDist1 = binCenters(minIdx1);
closestDist2 = binCenters(minIdx2);

% Now find all the pixels in pxls that correspond to these distances
tolerance = 0.01;       

% Find indices where the distance is close to the identified peaks within the tolerance
color1Indices = abs(dis_left - closestDist1) <= tolerance;
color2Indices = abs(dis_left - closestDist2) <= tolerance;

% Use logical indexing to get back the RGB values from the original pxls array
% that correspond to the remaining pixels after thresholding
color1RGB = remainingPxls(color1Indices, :);
color2RGB = remainingPxls(color2Indices, :);
color1RGB_mean = mean(color1RGB, 1);
color2RGB_mean = mean(color2RGB, 1);

%% DISPLAY INITIAL SELECTION
% Display the first color
subplot(1, 3, 1); % This creates a subplot in a 1x2 grid, at the first position
imshow(repmat(reshape(color1RGB_mean, 1, 1, 3), 100, 100)); % Scale the square size as needed
title('Color 1');

% Display the second color
subplot(1, 3, 2); % This creates a subplot in a 1x2 grid, at the second position
imshow(repmat(reshape(color2RGB_mean, 1, 1, 3), 100, 100)); % Scale the square size as needed
title('Color 2');

subplot(1, 3, 3); % This creates a subplot in a 1x2 grid, at the second position
imshow(repmat(reshape(edgeVal, 1, 1, 3), 100, 100)); % Scale the square size as needed
title('Color 3');

% Parse out which is brown and which is blue
img= imgCurrent;
if min(color1RGB_mean(3)-color1RGB_mean(1:2)) < min(color2RGB_mean(3)-color2RGB_mean(1:2))
    blueValu = color1RGB_mean;
    brownValu = color2RGB_mean;
else
    blueValu = color2RGB_mean;
    brownValu = color1RGB_mean;
end
%% We now have a preliminary color basis. However, its likely terrible, so next we improve it
colorsRGB = rescale(double([edgeVal; blueValu; brownValu]), 0,1);

% a basis with good gamut coverage usually has high saturation, and to a certain extent
% high luminence as well. So we correct this here:
colors_hsv = rgb2hsv(colorsRGB);
colors_hsv(3,2) = max(0.9, colors_hsv(3,2)); % focus on the blue because 
colors_hsv(3,3) = max(0.9, colors_hsv(3,3));
% convert back to RGB
saturatedColorsRGB = hsv2rgb(colors_hsv);
saturatedColorsRGB = min(max(saturatedColorsRGB, 0), 1);  % Clamp values to [0, 1]

%% CHROMATIC ADAPTATION
% Now we convert to the XYZ space, which is the appropriate space for the chromatic adaptation 
whiteVal_xyz = rgb2xyz(saturatedColorsRGB(1, :));
brownVal_xyz = rgb2xyz(saturatedColorsRGB(2, :));
blueVal_xyz = rgb2xyz(saturatedColorsRGB(3, :));

    % normalize white
whiteVal_xyz = whiteVal_xyz./sum(whiteVal_xyz(:));
% h = sort(whiteVal_xyz);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
% midWhite  = h(2);

%% estimate the equal energy estimate for new primaries
eqEnergyPrimaryEst = max([0 0 0], [1 1 1] - (brownVal_xyz + blueVal_xyz));

a1 = find(blueVal_xyz==max(blueVal_xyz));
a2 = find(brownVal_xyz==max(brownVal_xyz));
maxInds = setdiff([1 2 3], [a1 a2]); 
if length(maxInds)==2
    maxInds = maxInds(1);
end
a3 = find(eqEnergyPrimaryEst==max(eqEnergyPrimaryEst));
if ~ismember(a3, maxInds)
    k = max([eqEnergyPrimaryEst(a3), blueVal_xyz(maxInds), brownVal_xyz(maxInds)]);
    k2 = (1-k)/2+k;
    eqEnergyPrimaryEst(maxInds) = k2;
end
if eqEnergyPrimaryEst(maxInds)<0.49
    eqEnergyPrimaryEst(maxInds) = (1+eqEnergyPrimaryEst(maxInds))/2;
end

%% CALCULATE matrix M (still mostly following Bruce Lindbloom
partM = [brownVal_xyz; blueVal_xyz; eqEnergyPrimaryEst]';
B = whiteVal_xyz';  % Assuming midWhite normalization is not necessary

% If partM is not full rank, use pseudo-inverse
if rank(partM) < 3
    C = pinv(partM) * B;
else
    C = linsolve(partM, B);
end

% Construct matrix M
M = partM * diag(C);

%% apply transformation matrix M to the image!

[blue, brown, third] = applyColorMatrixTransformation(imgCurrent, M);

mixedData = cat(3, third, brown, blue);
dataForPy = py.numpy.array(mixedData);
ica = py.sklearn.decomposition.FastICA(pyargs('n_components', 3, 'random_state', 0));
S_ = ica.fit_transform(dataForPy);
independentComponents = double(S_);
mixingMatrix = double(ica.mixing_);

recoveredData = mixedData / mixingMatrix';

% Reshape the recovered data back into an image format if needed
recoveredImages = reshape(recoveredData, [numRows, numCols, numChans]);

figure;
imshow(vertcat(horzcat(imgCurrent(:,:,1),imgCurrent(:,:,2),imgCurrent(:,:,3)), horzcat(third, brown, blue)));

figure;
figure; montage({rescale(imsubtract(brown, blue),0,1), rescale(imsubtract(blue, third),0,1), rescale(imsubtract(third, brown),0,1)});
end



function thresh = findValley(matrix)

mat = matrix(:);
[counts, edges] = histcounts(mat, 'BinMethod','fd'); % this BinMethod is well suited for highly skewed data
dCounts = diff(counts); % derivative of counts

% finds the ties when deriv changes from + to -
valleyIndices = find(dCounts(1:end-1) < 0 & dCounts(2:end) >= 0, 1);

if isempty(valleyIndices)
    error('No valley found');
end

valleyIndex = valleyIndices(1);
thresh = edges(valleyIndex);
end





function [blue, brown, thirdChnl] = applyColorMatrixTransformation(img, M)
imXYZ = rgb2xyz(img);

pixXYZ = reshape(imXYZ, [],1,3);
pixXYZ = squeeze(pixXYZ).*100;
pixXYZ = pixXYZ';

% if the matrix has huge condition #, inverse fcn will not be accurate so we use the
% pseudoInverse
if cond(M) < 1e6
    newColorSpace = M \ pixXYZ;
else
    p = pinv(M); % less computationally demanding and more accurate for an  illconditioned mat
    newColorSpace = p * pixXYZ;
end


newSmall = reshape(newColorSpace', [size(img, 1:2), 3]);

blue = mat2gray(newSmall(:,:,1));
brown = mat2gray(newSmall(:,:,2)); 
thirdChnl = mat2gray(newSmall(:,:,3));
end


