function finMask  = vertexTweaking_handdrawn(outOfBoundsMask,imgCurrent, varargin)
% I made this funciton up as a last ditch effort to refine a
% segmentation effort. The mask is converted into a more simple shape with
% fewer verticies using the reducepoly function, the sensativity for which
% we carefull calibrate so that there is never more than some arbitrary number 
% of nodes, e.g. 125 points. User then can drag each vertex around to alter segmentation. 

% hidden behind the GUI is ANOTHER GUI with a prompt asking if you
% finished yet. After that, to refine further, a small dilating and a
% couple iterations of active contour are initiated, s.t it doesn't
% look like jagged and handdrawn. 

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
            tol = tol*1.5;
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
        message= {'When you are happy with mask, press OK'};
        uialert(alertFig, message, 'Confirm',...
            'Icon', 'question', 'Modal', false,...
            'CloseFcn', @buttsCallback);
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

function buttsCallback(~, ~)
    
    uiresume;
end