function finalMask = LinkUpBrokenEdges(edgeImage)
% finalMask = LinkUpBrokenEdges(edgeImage) INPUT MUST BE BW

% Now the strategy is to connect an endpoint to the closest other endpoint.
% There are two options.
%   (1) One option is whether or not to connect to other endpoints ONLY on a DIFFERENT curve, or,
%       if connecting to the other endpoint on the same curve is OK.
%   (2) The other option is how close the other endpoint must be.  It is a user specified distance in pixels.
%       It should be infinity if it will jump gaps regardless of how far away they are.
% Also in the strategy is that if an endpoint has been "used" so that it's no longer an endpoint,
% it should be excluded from connecting to any other end point.  We don't want two, three, or more lines emanating from an original endpoint.
% So we need to construct a list where every row has two endpoints, the labels (ID #'s) of the curves each endpoint is on,
% the distance between them, and whether the endpoints have been used (connected) or not.
% So we'll construct a table with these columns
% (x1, y1)    (x2, y2)    point1#    point2#    label1    label2    distance    used

%%
longestGapToClose = 80;

%%

% Find the endpoints
endPs = bwmorph(edgeImage, 'endpoints');
[endP_Rows, endPointCols] = find(endPs);
numberOfEndpoints = length(endP_Rows);

% Label the image.  Gives each separate segment a unique ID label number.
[labeledImage, ~] = bwlabel(edgeImage);

% Get the label numbers (segment numbers) of every endpoint,
% because we don't want to link an endpoint to the other endpoint on the same line.
theLabels = zeros(numberOfEndpoints, 1);
for k = 1 : numberOfEndpoints
    thisRow = endP_Rows(k);
    thisColumn = endPointCols(k);
    % Get the label number of this segment
    theLabels(k) = labeledImage(thisRow, thisColumn);
end

finalMask = edgeImage;	% Create/Initialize output image.  
lineCounter = 1;	% Initialize a counter of all the lines we draw no matter where they start and end.

f1 = uifigure('Visible', 'off'); 
ax1 = axes(f1);
imshow(edgeImage, []);

for k = 1 : numberOfEndpoints
    thisRow = endP_Rows(k);
    thisColumn = endPointCols(k);
    
    % Get the label number of this segment
    thisLabel = theLabels(k);
    
    % Get indexes of the other end points.
    otherEndpointIndexes = setdiff(1:numberOfEndpoints, k);
    
    otherLabels = theLabels(otherEndpointIndexes);
    onSameSegment = (otherLabels == thisLabel); % List of what segments are the same as this segment
    otherEndpointIndexes(onSameSegment) = []; % Remove if on the same segment
    
    % Now get a list of only those end points that are on a different segment.
    otherCols = endPointCols(otherEndpointIndexes);
    otherRows = endP_Rows(otherEndpointIndexes);
    
    % Compute distances
    distances = sqrt((thisColumn - otherCols).^2 + (thisRow - otherRows).^2);
   
    % Find the min - the one closest point.
    [minDistance, indexOfMin] = min(distances);
    nearestX = otherCols(indexOfMin);
    nearestY = otherRows(indexOfMin);
    
    if minDistance <= longestGapToClose
    
        % Draw line from this endpoint to the other endpoint.
        hLine = images.roi.Line(ax1, 'Position',[thisColumn, thisRow; nearestX, nearestY], 'LineWidth', 5);
        mask = createMask(hLine, edgeImage);

        
        finalMask = finalMask | mask;	% Burn line into output image.
        
        cc = bwconncomp(finalMask);
        if cc.NumObjects ==1
            continue
        end
       
        % Increment the line counter.
        lineCounter = lineCounter + 1;
    end
    
end
% f1.Visible = 'on';

% figure;
% imshowpair(edgeImage, finalMask);

%CHECK
% bwcon = bwconncomp(finalMask, 8);
% if bwcon.NumObjects ~=1
%     disp(strcat('there are : ', num2str(bwcon.NumObjects), ' in you new mask?'));
% end


