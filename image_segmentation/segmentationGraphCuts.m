function BWfull = segmentationGraphCuts(imgCurrent)
% Segment an RGB color image into foreground and background, using
% iterative graphi
% NOTE: for some optimization of this function, you'll need to change some
% parameters located at the top of this code. Anyways, this function
% covers the preprocessing, and postprocessing of the builtin function
% GRABCUTS. 

% Michael Glendinning, 2022

img = im2double(imgCurrent);
img = padarray(img, [1,1,1], 1, 'both');

img(img(:,:,:)>0.9)=1;
img = img.^4; %will darken everything except the 1's which don't change. values close to 1 get less dark proportionally

[M,N] = size(img);

%% specification of some parameters

howManyPartitionsPerShortSide = 16; %this variable is a means of specifying how
% densely you want points to be added to the poly line you draw, 
% intersdisperesed evenly between your nodes. The more points given to
% graphCuts, the better. 


sizeOfSuperPixel = 100; %this variable identifies the extent to which you want
% to oversegment the input image to, a requirement of the grab cuts
% algorithm.  100 refers to the approximate size of each superpixel area. 
% specifically: sizeOfSuperPixel = prod(size(myImg,1:2))/(# desired superpixels)

%% start processing
]]
imgShrunk = imresize(imgCurrent, 0.5, 'lanczos3');
distBetweenPoints = min(size(imgShrunk, 1:2))/howManyPartitionsPerShortSide;

if ndims(imgShrunk)<3
    imgShrunk = cat(3, imcomplement(imgShrunk), zeros(size(imgShrunk, 1:2), 'double'), zeros(size(imgShrunk, 1:2), 'double'));
end
    
f1 = uifigure;
ax1 = axes('Parent', f1);
drawnow;
imshow(imgShrunk, 'Parent', ax1);

%first the foregrounds
disp('the first line to draw is FOREGROUND. The closer the line is to the edge, the better');
disp('double-click to finish');
inboundsPre = drawpolyline(ax1, 'Color', [0.4660 0.6740 0.1880]);

wait(inboundsPre);
set(f1, 'Visible', 'off');
newVerticesIn = fillInPolyLine(inboundsPre,6, distBetweenPoints);

%graph based segmentation works better without lots of redundant pixel
%data.... 
P = impixel(imgShrunk, newVerticesIn(:,2), newVerticesIn(:,1));
[~, indx] = unique(round(P, 3), 'rows');
myrange = 1:1:length(P);
newVerticesIn(~ismember(myrange, indx), :) = [];

delete(inboundsPre); pause(0.1);
inbounds = drawpolyline(ax1, 'Color', 'g', 'Position', newVerticesIn, 'LineWidth', 1);
set(f1, 'Visible', 'on');
inBoundsPos = wait_adjustPolyline_getFinalCoords(inbounds);
delete(inbounds);

%next is the background
outOfBoundsPre = drawpolyline(ax1, 'Color', [0.6350 0.0780 0.1840]);
newVerticesOut = fillInPolyLine(outOfBoundsPre,6, distBetweenPoints);

%as above
P2 = impixel(imgShrunk, newVerticesOut(:,2), newVerticesOut(:,1));
[~, indx2] = unique(round(P2, 2), 'rows');
myrange2 = 1:1:length(P2);
newVerticesOut(~ismember(myrange2, indx2), :) = [];

delete(outOfBoundsPre); pause(0.1);
outOfBounds = drawpolyline(ax1, 'Color', 'r', 'Position', newVerticesOut, 'LineWidth', 1);
outBoundsPos = wait_adjustPolyline_getFinalCoords(outOfBounds);
delete(outOfBounds);
close all force
%% 
inBoundsInd = round(sub2ind(size(imgShrunk),inBoundsPos(:,2),inBoundsPos(:,1)));
OutBoundsInd = round(sub2ind(size(imgShrunk),outBoundsPos(:,2),outBoundsPos(:,1)));

numSuperPix = floor(prod(size(imgShrunk,1:2))/sizeOfSuperPixel);

f1 = figure('Visible', 'on'); 
imshow(imgShrunk);
drawnow;

x = inBoundsPos(:,1); y = inBoundsPos(:,2);
k = boundary(x, y, 0.25);

mypolygon = drawpolygon('Position', [x(k),y(k)]);
myMask = imdilate(createMask(mypolygon), strel(ones(180))); %
hold on
visboundaries(myMask);
Lrgb = superpixels(imgShrunk, floor(numSuperPix/2));
BW = grabcut(imgShrunk, Lrgb, imdilate(createMask(mypolygon), strel(ones(160))),inBoundsInd,OutBoundsInd);
visboundaries(BW); 
hold off

BW2 = imclose(BW, strel('octagon', 7));
BW3 = imdilate(imfill(BW2, 'holes'), strel('disk', 15, 8));
BW4 = activecontour(rgb2gray(imgBW), BW3, 25, 'edge','ContractionBias', 1);

BWfull = imresize(BW4, size(imgCurrent, 1:2));

end


function coordMatrix = wait_adjustPolyline_getFinalCoords(hROI)
% this function is sort of shit way to accomplish this but oh well

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);
uiwait; % sister uiresume is hidden in listener
delete(l); % Remove listener

% Return the current coordinates of points of polyLine
% coordMatrix will be n-by-2 2d matrix, ordered list of rows of 2d coords 
coordMatrix = hROI.Position;
end

function clickCallback(~,evt)
% this listener prevents anything from moving fwd until you double click
if strcmp(evt.SelectionType,'double') 
    uiresume;
end
end   

