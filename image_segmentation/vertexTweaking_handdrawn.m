function finMask  = vertexTweaking_handdrawn(outOfBoundsMask,imgCurrent, varargin)
% syntax: finMask  = vertexTweaking_handdrawn(outOfBoundsMask,imgCurrent) OR
%         finMask  = vertexTweaking_handdrawn(outOfBoundsMask,imgCurrent, NUMBER_OF_NODES)
% When automated refinements fall short, when permissible at least, one can turn to manual
% refinement. This funciton converts a binary mask into a polygon shape (only 1). 
% # verticies is limited using the reducepoly function, the sensativity for which
% we can calibrate with a while-loop to match a requested input number of nodes.
% If not supplied, function assumes value of 125 points. After polygon is drawn, 
% it is overlaid on the RGB image, and then the user then can drag each vertex 
% around to alter segmentation. A second GUI should also be visible giving you the option
% to indicate you're finished refining segmentation and would like to submit. 

% Just drag and drop vertices to move. After submission, to smooth any jagged edges 
% created as a result of polygon simplification, a small dilation is applied, followed
% by a couple iterations of active contour (i.e. snakes). 

% note, number of nodes is approximate, may vary by as much as 10

% -Michael Glendinning, 2023

if isnumeric(varargin{1})
    approx_nPoints = varargin{1};
else
    approx_nPoints = 125;
end

   close all force
        outOfBoundsMask = bwareafilt(outOfBoundsMask,1);
        boundaries = bwboundaries(outOfBoundsMask,8); 
        largest = boundaries{1};
        largestSize = size(boundaries{1},1);
     
        for kk = 2:numel(boundaries)
            newSize = size(boundaries{kk}, 1);
            if newSize>largestSize
                largestSize = newSize;
                largest = boundaries{kk};
            end
        end
        % we start with a very low tolerance (ie huge # of points) and
        % gradually get more restrictive until we have like 100 points
        % (#s arbitrarily chosen for my data)
        tol = 0.0005;
        vert2 = reducepoly(largest, tol);
        while size(vert2, 1)>approx_nPoints
            tol = tol*1.4; %arbitrary scale factor (shrinking this value closer to 1 will increase precision, at the cost of efficiency)
            vert2 = reducepoly(largest, tol);
        end
        x = vert2(:,2);
        y = vert2(:,1);
        vert3 = [x, y];
        
        f2 = figure;
        movegui(f2, 'east');
        ax2 = axes(f2);
        imshow(imgCurrent,'Parent',ax2, 'Border', 'tight');
        hold on
        roi = drawpolygon('Position', vert3, 'Color', 'r');
        hold off
        
        alertFig = uifigure;
        movegui(alertFig, 'west');
        message= {'When you are happy with mask refinements, press OK'};
        uialert(alertFig, message, 'Confirm',...
            'Icon', 'question', 'Modal', false,...
            'CloseFcn', @buttCallback);
        uiwait;
        
        outOfBoundsMask_temp = createMask(roi, imgCurrent);
        drawnow;
        
        close all force
        
        se3 = strel('disk', 7, 8);
        outOfBoundsMaskbig = imdilate(outOfBoundsMask_temp, se3);
        
        if ndims(imgCurrent) == 3
            imgCurrent = imadjust(rgb2gray(imgCurrent));
        end
        
        finMask = activecontour(imgCurrent, outOfBoundsMaskbig, 16, 'edge');

end

function buttCallback(~, ~)
    uiresume;
end
