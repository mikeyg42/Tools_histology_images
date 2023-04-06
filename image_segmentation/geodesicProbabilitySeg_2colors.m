
function outOfBoundsMask = geodesicProbabilitySeg_2colors(imgCurrent)
% syntax outOfBoundsMask = geodesicProbabilitySeg_2colors(imgCurrent)
    % This calls the IPT function for the geodesic soft-segmentation method
    % (a "distance-based color segmentation algorithm"), optimizing for WSI and 
    % semi-automating the requisite inputs. 
    
    % Notably, the builtin IPT function I've found outputs a slightly dilated
    % mask... I usually have better segmentation using the "probability map" that
    % is an optional auxillary output. I threshold that map at 95% level,
    % and then multiply it with the normal output mask to mitigate any dilation. 

    % Michael Glendinning, 2022
    
 close all force
    
% 1st step : make image smoother, as homogenous areas segment easer. 
% I'm opting here for the edge-preserving imbilateral filter
    
% define a small central patch to determine im variance
  sz = size(imgCurrent, 1:2);
  topLeftCorner = sz./2-50;
  rect = [topLeftCorner(2), topLeftCorner(1), 100, 100];
  patch = imcrop(imgCurrent, rect);
        
% these 3 lines are directly from the documentation for IMBILATFILT
%      ||
%      vv 
 edist = patch.^2;
 edist = sqrt(sum(edist,2)); % Euclidean distances from origin
 patchVar = var(edist(:));
        
 currentImageAdj = imbilatfilt(imgCurrent, patchVar*3, 9); %this fcn is rather slow
        
% 2nd: the target foreground regions is indicated by user via their making of  "scribbles". 
% the geodesic segmentation builtin fcn takes global properties of the defined 
% region, using built-in gabor filtures, s.t. the larger the region the better! 
% BUT, it won't work if you cross your ROIs or borders. 
  
   % I've already pre-defined the background region, so only the foreground
   % needs manual "scribbles" (ie a polygon >100 pixles )
   
RI = imref2d(size(currentImageAdj,1:2));
f1 = figure('Visible', 'off'); 
ax1 = axes(f1);
imshow(currentImageAdj, RI, 'Parent', ax1);
        
   % this draws a predefined "picture frame" roi circumscribing the image extents, is 12 wide
polyPos = [1 1; size(currentImageAdj,2) 1; size(currentImageAdj,2) size(currentImageAdj,1);...
           1 size(currentImageAdj,1); 1 12; 12 12;12, (size(currentImageAdj,1)-12);...
           (size(currentImageAdj,2)-30) (size(currentImageAdj,1)-30);...
           (size(currentImageAdj,2)-12), 12;11 10; 1 1];
        
%set(f1, 'Visible', 'on')

disp(' draw a polygon around a large area representative of the foreground');
disp('(i.e. free of debris, rips, or anything like that.) Note: polygon does NOT need to span the tissue...'); 
disp(' your goal is to select a representative sample... when done double-click on a node to submit');

hold on
%My background region is immutable (as far as the user is concerned)
hBorder = drawpolygon('Parent', ax1, 'Position', polyPos, 'LineWidth', 2, 'InteractionsAllowed', 'none', 'Color', 'y');

%allows user to draw scribbles covering extent of foreground region
hForeground = drawpolygon('Parent', ax1, 'LineWidth', 2, 'Color', 'm', 'InteractionsAllowed', 'all'); 
wait(hForeground);
        
%NOTE : double-clicking on the ROI indicates that you are done interacting with the ROI.   
hold off
set(f1, 'Visible', 'off');
        
maskBackground = logical(createMask(hBorder));
maskObject = logical(createMask(hForeground));
        
delete(hBorder);
delete(hForeground);
close all force;
        
% call the IPT geodesic segmentation function, with adaptive channel weighting
[labeledSegmentation,probabilityMap] = imseggeodesic(currentImageAdj,maskObject,maskBackground, 'AdaptiveChannelWeighting', true);
binarySeg = labeledSegmentation;
val =labeledSegmentation(2,2);

if val == 2
      binarySeg(binarySeg==2)  = 0;
else
      binarySeg(binarySeg== 1) = 0;
      binarySeg = binarySeg./2;
end

mask = logical(binarySeg);
            
% special case for when k = 2 colors - the probability maps are inverse images of each other.
[pMap, ~] = imsplit(probabilityMap);
      
pMap(pMap>0.95) = true; %select anything that has more than 95% confidence
pMap(pMap<0.95) = false; 
        
mask = pMap.*mask;
         
% fill in any holes! (optional, mainly for rips/tears in tissue, ect.)
outOfBoundsMask = imfill(mask, 8, 'holes');
          
end
