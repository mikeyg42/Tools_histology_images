function newVertices = fillInPolyLine(original_polyline, MIN_distBewtweenPoints, MAX_distBetweenPoints)
% returns a new vector of points defining a line in which all points are within a certain
% specified distance of each other, defined as MAX_distBetweenPoints. 
% All edges longer than MAX_distBetweenPoints threshold are split up.

% original_polyline can be either a polyline ROI, or it can be a matrix of
% vertices

% if no minimum is desire, just put 0.
% if no maximum is desired, just put Inf

% if you want to later create a new ROI using new points, call: 
%      newPolyLine = drawpolyline('Position',newVertices);

% -Michael Glendinning

try MyVertices = original_polyline.Position; % verticies in the form [x, y]
    
catch ME
    if ismatrix(original_polyline) && size(original_polyline, 2)==2
        MyVertices = original_polyline;
    else 
        rethrow(ME);
    end
end

%%
format longG

% to avoid index issues I split this up into three seperate loops

%% loop 1 removes colinear points! ################################

% I should change this to calculate the determinant of this matrix 
%[x1, y1, 1; x2,y2,1;x3,y3,1] which is 0 if colinear
numPoints = size(MyVertices,1);
connectionMatrix = zeros(numPoints,1);

    for ii=2:numPoints
        v1 = MyVertices(ii-1,:);
        v2 = MyVertices(ii,:);
        
        if any((v2-v1)==0)
        connectionMatrix(ii, 1) = connectionMatrix(ii-1, 1)+1;
        else
           connectionMatrix(ii-1, 1) = 0 ;
       end
    end

MyVertices = MyVertices.*(connectionMatrix==0);
MyVertices(all(~MyVertices,2), : ) = [];

%% loop two address the minimum distance between points constraint #########
if MIN_distBewtweenPoints > 0
numPoints = size(MyVertices,1);
edgeLength = zeros(numPoints,1);

    for ii=2:numPoints
        v1 = MyVertices(ii-1,:);
        v2 = MyVertices(ii,:);
        edgeLength(ii,1) = pdist([v1;v2],'euclidean');
    end
    edgesToDelete = find(edgeLength(:,1)<MIN_distBewtweenPoints); 
    edgesToDelete(1) = []; % 1 will always be on the list because it will have edgelength = 0 bc that row is never populated by the loop! 
    vertsToDelete = zeros(length(edgesToDelete),1);
    % edges have 2 vertices, so i need to decide which of the 2 to delete.
    % I've opted to prioritize a smaller boundary.
    
    % ...We don't want to lose first or second points, 

    for jj = 1:length(edgesToDelete)
        if edgesToDelete(jj, 1) < 3 
            vertsToDelete(jj,1) = 2; % and not 1, I never want to delete endpoints!
            continue
        else
            vert0 = MyVertices((edgesToDelete(jj, 1)-2),:);
        end
        
        vert1 = MyVertices((edgesToDelete(jj, 1)-1),:);
        vert2 = MyVertices(edgesToDelete(jj, 1),:);
        
        if edgesToDelete(jj, 1)+1>length(edgeLength) 
            vertsToDelete(jj,1) = length(edgeLength)-1; %again preserving end points!
            continue
        else
        vert3 = MyVertices(edgesToDelete(jj, 1)+1,:);
        end
        
        norm1 = pdist([vert0;vert2],'euclidean');
        norm2 = pdist([vert1;vert3],'euclidean');
        
        if norm1>norm2
        vertsToDelete(jj,1) = edgesToDelete(jj, 1);
        elseif norm1<norm2 || norm1==norm2
        vertsToDelete(jj,1) = edgesToDelete(jj, 1)-1;
        end
    end
    vertsToDelete = unique(vertsToDelete, 'rows');
    for kk = 2:1:length(vertsToDelete)-1
        if vertsToDelete(kk)+1==vertsToDelete(kk+1) && vertsToDelete(kk)-1==vertsToDelete(kk-1)
        vertsToDelete(kk) = [];
        end
    end
    
    MyVertices(vertsToDelete,:)=[];
end



%% loop number 3 determines how many points to add in  ##################
numPoints = size(MyVertices,1);
totalNewPoints = 0;
numberNewPointsADDED = zeros(numPoints,1);

for ii=2:numPoints
v1 = MyVertices(ii-1,:);
v2 = MyVertices(ii,:);

% equivilent to dist = norm(v1-v2)  ;
dist = pdist([v1;v2],'euclidean');

numberNewPointsADDED(ii,1) = floor(dist/MAX_distBetweenPoints);
totalNewPoints = totalNewPoints+numberNewPointsADDED(ii,1);
end

%% PSYCH! loop #4... assigns the points chunk by chunk into output format ##########

%preallocate new vertex matrix
totalPoints = totalNewPoints + numPoints;
newVertices = zeros(totalPoints,2);

%prepopulate the first row with first vertex
newVertices(1,1:2) = MyVertices(1,1:2);

lastFullRowInd = 1;
%%
for jj=2:numPoints

x1 = MyVertices(jj-1,1);  y1 = MyVertices(jj-1,2);
x2 = MyVertices(jj,1);  y2 = MyVertices(jj,2);

linspaceLength = numberNewPointsADDED(jj,1)+2; %linspace includes your start and end as points
yPoints = linspace(y1, y2, linspaceLength);
xPoints = linspace(x1, x2, linspaceLength);
newMatrixPointsFill = [xPoints', yPoints'];

%remove the top row because it was already assigned in previous iteration of loop
newMatrixPointsFill(1, :) = [];

startFill = lastFullRowInd+1;
endFill = startFill+linspaceLength-2;
newVertices(startFill:endFill, :) = newMatrixPointsFill;

lastFullRowInd = endFill;

end
newVertices = unique(newVertices, 'rows');
end

