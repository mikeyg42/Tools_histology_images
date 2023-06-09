function  choosingROIs
close all force
warning('off', 'all');

%% the preamble

% call the function into which you've inputted all the user preferences and datasets to be
% here used.
mySettings = setts_and_prefs;

% call your parseDataset function, which interprets settings and loads up the dataset into
% a structure called "data"
data = parseDataset(mySettings, 'choosingROIs'); % inside data will be a image datastore from which you can read in images

% load up and check over your input files
imgCurrent = readimage(data.rgbIMDS, 1);
imgCurrent = im2double(imgCurrent);
maskCurrent = readimage(data.maskIMDS, 1);
if size(imgCurrent, 1:2)~= size(maskCurrent)
    maskCurrent = imresize(maskCurrent, size(imgCurrent,1:2), {@oscResampling, 4});
end
totalBlobArea = sum(imcomplement(maskCurrent(:)));

%% start of the processing to generate 10 ROI's
close all

%initialize mask for storing locations and counter of ROIs
numSubsets=1; % the counter fr while loop
FINmasks{mySettings.chooseROI.numROIs} = [];


% Calls the initialization routine, which comes up with a good 1st pass at a solution
[solution1, tracking, info] = circlePacking_initiate(maskCurrent, mySettings.chooseROI.sizeROI, 3, mySettings);
solution1(1).all = solution1(1).all | ~maskCurrent ;

% save this result into a cell array you'll populate each iteration of the main loop, for evaluation of processes at the very end
finishedLoops{1} = solution1(1).all;

nCircles = nnz(unique(round(tracking(:,1), 2)));

% these 2 arrays always we actively update so it includes whats CURRENTLY in the mask.
currentRadii  = tracking(:, 2);
currentRadii = currentRadii(currentRadii ~= 0);
currentCenters(:, 1:2) = tracking(:, 3:4);
currentCenters = reshape(currentCenters(currentCenters ~= 0), [], 2);

% Construct unweighted, undirected graph
circlegraph = makeCircleNetworkGraph(currentRadii, currentCenters, maskCurrent);

% identify those nodes which have only 1 or 2 alters
[id, alter1, alter2] = select3Circles2reposition(circlegraph, currentCenters);

% initiallize
repositionedCenter = [];
counting = 1;

%% Start Loop - note, first pass through the loop skips this first if-then block which recapitulates the code above mostly
while counting<6

    % repositionedCenter ONLY == 1 in the first loop. Always after that it will be > 1
    if counting > 1
        % recreate these clean using the "winning"arrays
        currentCenters = winningCenters{counting-1};
        currentRadii = winningRadii{counting-1};

        nCircles = numel(currentRadii); % needs to be re-evaluated with each loop!

        % recalculate the undirected graph of connectivity
        circlegraph = makeCircleNetworkGraph(currentRadii, currentCenters, maskCurrent);

        [id, alter1, alter2] = select3Circles2reposition(circlegraph, currentCenters);
    end
   
    % Check that the binary image you are about to fill in the holes of didnt get inverted
    unfilledImage = ~finishedLoops{counting};
   test = imclearborder(unfilledImage);
   if ~isequal(~finishedLoops{counting}, test)
       unfilledImage = imcomplement(unfilledImage);
   end

    % make the mask with 3 circles filled back in.
    mask1 = solution1(counting).individ{alter1} ;
    mask2 = solution1(counting).individ{alter2} ;
    mask3 = solution1(counting).individ{id} ;
    filledInHoles = (unfilledImage | mask1 | mask2 | mask3);
   drawnow;

    if counting == 1 % I've yet to figure out why this is needed but for now... it stays or the code breaks!
        filledInHoles = imcomplement(filledInHoles);
    end

    %  Add these points to the repositioned list
    [y1,x1] = find(mask1);[y2,x2] = find(mask2);[y3,x3] = find(mask3);
    cPoints = [mean(y1(:)), mean(x1(:)); mean(y2(:)), mean(x2(:));mean(y3(:)), mean(x3(:))];
    [~, iIi] = pdist2(currentCenters, cPoints, 'Euclidean', 'Smallest',3);
    iii = unique(reshape(iIi', [], 1), 'rows', 'stable');

    repositionedCenter = [repositionedCenter; currentCenters(iii(1:3), :)];

    % apply the covering up mask to the previous solution to create NEWCHALLENEGE mask
    newChallenge = ~filledInHoles; % white background, back unfilled area, white circles ontop
    if isfloat(newChallenge)
        newChallenge = logical(newChallenge);
    end

    % sanity check: confirm that we didn't get inverted inadvertently
    test = imclearborder(newChallenge);
    if ~isequal(test, newChallenge)
        newChallenge = imcomplement(newChallenge);
    end
    % the following COMPLETELY relies on newChallenge NOT being inverted (it needs black
    % background and white foreground

    %use small morphological filtering to ensure no hangnails or extra bits tag along
    workingBlob = newChallenge & imdilate(bwareafilt(imerode(newChallenge, strel(ones(7))),1), strel(ones(7)));

    % This script is the actual packing repacking script. We call both variations 1-by-1 and then pick the best!
    [masks1, cents1, rads1] = fillInBlob(workingBlob, info, maskCurrent, 1); drawnow;

   masks1mat = reshape(cell2mat(masks1), size(newChallenge,1),size(newChallenge,2),[]);
   maskPacked_type1 = newChallenge & ~logical(sum(masks1mat, 3));

   deltaTheta = queryBlobs(masks1, cents1, rads1, maskPacked_type1, maskCurrent, info); 



    [masks2, cents2, rads2] = fillInBlob(workingBlob, info, maskCurrent, 2); drawnow;
     masks2mat = reshape(cell2mat(masks2), size(newChallenge,1),size(newChallenge,2),[]);
    maskPacked_type2 = newChallenge & ~logical(sum(masks2mat, 3));



    mergedMasks = {maskPacked_type1, maskPacked_type2};

    val1 = sum(maskPacked_type1(:));
    val2 = sum(maskPacked_type2(:));
    [arealeft, ix_val] = min([val1, val2]);

    if ix_val == 1
        disp( 'Chosen direction of packling is TOP DOWN');
    else
        disp( 'Chosen direction of packling is BOTTOM UP');
    end

    cents_all = {cents1, cents2};
    rads_all = {rads1, rads2};

    currentCenters = [currentCenters; cents_all{ix_val}];
    currentRadii = [currentRadii; rads_all{ix_val}];

    masks_all = {masks1, masks2};
    ALLMASKS = [solution1(counting).individ, masks_all{ix_val}];

    winner = mergedMasks{ix_val};

    isoBlobs = imclearborder(imcomplement(imdilate(winner, strel('disk', 37, 8))));
    blobProps = regionprops('table',isoBlobs, 'Centroid');
    distError = pdist2(blobProps.Centroid, currentCenters, 'Euclidean');
    [~, idx] = min(distError, [],2);

    winningCenters{counting} = currentCenters(idx, :);
    winningRadii{counting} = currentRadii(idx);

    solution1(counting+1).individ = ALLMASKS(idx);

    clear current*

    disp(strcat('Total packed : ' , num2str(numel(winningRadii)), ' with coverage of : ', num2str(1-arealeft/totalBlobArea)));

    newChallenge = winner;
    figure; imshowpair(isoBlobs, newChallenge, 'montage');

    counting = counting+1;
    % note: its not an error that this comes after updating "counting"... finishedLoops{1} is already full with the initialization fcn output!!
    finishedLoops{counting} = newChallenge;

end

figure;
montage(finishedLoops)


end

% ==========================================================================
% ==========================================================================

function graphOut = makeCircleNetworkGraph(currentRadii, currentCenters, tissueMask)
% graphOut = makeCircleNetworkGraph(currentRadii, currentCenters, tissueMask)

warning('off', 'all');

wiggleRoom = 0.1;

nCircles = size(currentRadii, 1);

% Use the sum of 2 circles' radii to define what the theoretical max distance 2 circles
% could be seperated (center to center) while still "kissing",as compared to real distances
[radius1, radius2] = meshgrid(currentRadii);
minSpanDist = radius1 + radius2;
minSeperation = minSpanDist * (1+wiggleRoom); % give 0.1% percentage of leniency and avoid rounding issue

actualSeperation = pdist2(currentCenters, currentCenters, 'euclidean');
withinMinSeparation = actualSeperation <= minSeperation;
withinMinSeparation = tril(withinMinSeparation, -1);
[Ith, Jth] = find(withinMinSeparation);

% Incorporate the outer Boundary as its own node
BWouter = bwboundaries(tissueMask, 8, "noholes",  'CoordinateOrder' , 'xy', 'TraceStyle', 'pixeledge');
outermaskxy = BWouter{1};

% Define distance outside of mask's boundary where the "node" will be
bufDist = 50; %don't make it so close that you start getting some hypotenuse issues. but not too far that you run into negative coord issues

outerPgon = polybuffer(polyshape(outermaskxy(:, 1), outermaskxy(:, 2), 'KeepCollinearPoints', true), bufDist, 'JointType', 'round');
outerDist = pdist2(outerPgon.Vertices, currentCenters, 'euclidean');
distances2out = min(outerDist, [], 1);

O_U_outside = distances2out' - bufDist - currentRadii(:);

outInd = O_U_outside < 10; % another way of providing some wiggly-wiggle room to wiggle
connect_outer = [repmat(nCircles + 1, nnz(outInd), 1), find(outInd)];

EdgeTable = table([Ith, Jth; connect_outer], 'VariableNames', {'EndNodes'});
graphOut = graph(EdgeTable);

if numnodes(graphOut) ~= nCircles+1
    graphOut = addnode(graphOut, nCircles+1-numnodes(graphOut));
end

end

% ==========================================================================
% ==========================================================================

function [masks, cents, rads] = fillInBlob(workingBlob, info, maskCurrent, topORbot)
% [masks, cents, rads] = fillInBlob(workingBlob, info, maskCurrent, UP_or_DOWN)

warning('off', 'all');
rads = [];
cents = [];
masks{1} = [];

counter = 1;
% with each loop of this while-loop, one circle is added, until we no longer can fit any
while counter < 11
    % regenerate the constrained triangulation given that we have modified our blob
    nHoles = 1 - bweuler(workingBlob);

    boundary2Fill = traceClump(workingBlob); %yx

    P = fliplr(boundary2Fill); %xy
    N = size(P, 1)-1; % last point restates the initial point, which we do not need
    C = [(1:(N-1))' (2:N)'; N 1]; %constraints the triangulation

    % If there are holes, we must include their boundaries for a proper triangulation
    if nHoles == 1
        % case of the unconnected circle...
        holePic = imclearborder(imcomplement(workingBlob));
        P = [P; fliplr(traceClump(holePic))];
        N2 = size(P, 1)-1; % last point restates the initial point, which we do not need
        C2 = [(N+1:(N2-1))' (N+2:N2)'; N2 N+1]; %constraints the triangulation
        C = [C; C2];
    elseif nHoles > 1
        disp('AHH i nOt sUrEE yEtTTT');
        return;
    end

    cdt = delaunayTriangulation(P, C);
    TF = isInterior(cdt);
    tr = triangulation(cdt.ConnectivityList(TF,:),cdt.Points(:,1),cdt.Points(:,2));

    %figure; imshow(workingBlob); hold on; triplot(tr);

    % A circumcenter is the center of the only circle one can draw connecting all 3 vertices of a triangle.
    [circums, radii] = circumcenter(tr);

    % define solution space by those circles with radius that fits in your range
    iC = radii<info.maxRadius & radii>info.minRadius;

    % EXIT CONDITION MAIN LOOP:
    % if there are no triangles left in ourn triangulation that fall within our search range,
    % then we have exhausted four search and need to leave!!
    if sum(iC)==0
        return
    end

    %define the solution space as all values pass through this size sieve
    solutionSpaceC = circums(iC, 1:2);
    solutionSpaceR = radii(iC, 1);

    % restrict solution space to those circles that the distance from center to the nearest
    % point on edge of blob is >= radius
    dcircum = pdist2(P, solutionSpaceC, 'Euclidean', 'Smallest',1);
    tf = dcircum' >= solutionSpaceR;
    solutionSpaceC(tf, :) = [];
    solutionSpaceR(tf) = [];
    if size(cents,1)>=1
        [restrict2, rI]  = pdist2(solutionSpaceC,cents, 'Euclidean', 'Smallest', 2*size(cents,1));
        toDestroy = sort(unique(rI(restrict2 < info.minRadius*2)), 'descend');
        if ~isempty(toDestroy)
            if numel(toDestroy) == size(solutionSpaceR, 1)
                return;
            end
            for h = 1:size(toDestroy, 1)
                solutionSpaceC(toDestroy(h), :) = [];
                solutionSpaceR(toDestroy(h)) = [];
            end
        end
    end

    %EXIT CONDITION #2 FOR MAIN LOOP
    if size(solutionSpaceR, 1) ==  0
        return
    end

    % Make absolutely sure that no points in the background/black will be selected
    boundaryAgain= traceClump(workingBlob);
    pgon2 = polyshape(boundaryAgain(:,2), boundaryAgain(:,1));

    TFin = isinterior(pgon2, solutionSpaceC(:, 1), solutionSpaceC(:, 2));
    solutionSpaceC(~TFin, :) = [];
    solutionSpaceR(~TFin, :) = [];

    % I've added here quick check of the nearest neighbor in the solution space to see if the neighbor is bigger...
    % Bigger should mean better -> a more connected packing! but not always. so we tread
    % carefully
    if size(solutionSpaceR,1)>2
        if topORbot ==1 %FILLS TOP DOWN
            [~ , nextPtIdx] = min(solutionSpaceC(:, 2)-solutionSpaceR(:), [], 'all');
            altNextPtIdx = knnsearch(solutionSpaceC, solutionSpaceC(nextPtIdx, :), 'k', 2);
            altIdx = altNextPtIdx(altNextPtIdx~=nextPtIdx);
            if  solutionSpaceR(altIdx)> solutionSpaceR(altIdx) && solutionSpaceR(altIdx)+solutionSpaceC(altIdx, 2)-solutionSpaceR(nextPtIdx)-solutionSpaceC(nextPtIdx, 2) > 0
                nextPtIdx = altIdx;
            end

        elseif topORbot ==2 %FILLS BOTTOM UP
            [~, nextPtIdx] = max(solutionSpaceC(:, 2)+solutionSpaceR(:), [], 'all');
            [d, altNextPtIdx] = pdist2(solutionSpaceC, solutionSpaceC(nextPtIdx, :),'Euclidean', 'Smallest', 2);
            altIdx = altNextPtIdx(2);

            if 2*(solutionSpaceR(altIdx)-solutionSpaceR(nextPtIdx))>d(2)
                nextPtIdx = altIdx;
            end

            try pgon3 = polybuffer(solutionSpaceC(nextPtIdx, :), 'points',solutionSpaceR(nextPtIdx));
                pgon = polybuffer(cents(end, :), 'points',rads(end));
                if overlaps(pgon3, pgon)
                    tf_in = isinterior(pgon, solutionSpaceC(:, 1), solutionSpaceC(:, 2));
                    solutionSpaceC(tf_in, :) = [];
                    solutionSpaceR(tf_in, :) = [];
                    [~ , nextPtIdx] = min(solutionSpaceC(:, 2)-solutionSpaceR(:), [], 'all');
                end
            catch
                [~ , nextPtIdx] = min(solutionSpaceC(:, 2)-solutionSpaceR(:), [], 'all');
            end
            drawnow;
        end
    elseif size(solutionSpaceR,1)==1
        nextPtIdx = 1;
    end

    if isempty(solutionSpaceR) || isempty(solutionSpaceC)
        return
    end

    if topORbot==1
        if size(solutionSpaceR,1)==2
            hROI1 = images.roi.Circle('Center', solutionSpaceC(1, :), 'Radius', info.maxRadius);
            hROI2 = images.roi.Circle('Center', solutionSpaceC(2, :), 'Radius', info.maxRadius);
            hROI3 = images.roi.Circle('Center', solutionSpaceC(1, :), 'Radius', solutionSpaceR(1,1));
            hROI4 = images.roi.Circle('Center', solutionSpaceC(2, :), 'Radius', solutionSpaceR(2,1));
        else
            radii2test = {info.maxRadius, (info.maxRadius + solutionSpaceR(nextPtIdx,1))/2, solutionSpaceR(nextPtIdx,1), info.minRadius};
            hROI1 = images.roi.Circle('Center', solutionSpaceC(nextPtIdx, :), 'Radius', radii2test{1});
            hROI2 = images.roi.Circle('Center', solutionSpaceC(nextPtIdx, :), 'Radius', radii2test{2});
            hROI3 = images.roi.Circle('Center', solutionSpaceC(nextPtIdx, :), 'Radius', radii2test{3});
            hROI4 = images.roi.Circle('Center', solutionSpaceC(nextPtIdx, :), 'Radius', radii2test{4});
        end
        ROIs = {hROI1, hROI2, hROI3, hROI4};
        % this loop assumes that bigger is better. So it starts off by checking out the largest
        % ROI (which is hROI1) and seeing if it meets the criteria of basically no overlap w/ workingBlob
        ROIindex = 1;

        while ROIindex<=4
            curROI = ROIs{ROIindex};
            tester = createMask(curROI, maskCurrent);
            mytest = tester.*workingBlob;

            drawnow;

            if sum(mytest(:))+25 >= sum(tester(:)) %this arbitrary 25 extra is just wiggle room to avoid spurrious rejections
                if size(solutionSpaceR,1)~=2
                    solutionSpaceR(nextPtIdx, :) = curROI.Radius;
                end

                [deltaCP, deltaRadius] = holeCheck(workingBlob, curROI.Center, curROI.Radius, maskCurrent);

                cents(counter,1:2) = curROI.Center+deltaCP;
                rads(counter, 1) = curROI.Radius+deltaRadius;
                masks{counter} = tester;

                workingBlob = bwareafilt(~tester & workingBlob, 1);
                ROIindex = Inf;

                delete curROI % this is probably unneceessary?
                drawnow;

            elseif sum(mytest(:))+25 < sum(tester(:)) && ROIindex == 4 % this means we exhausted all options and still we failed
                [deltaCP, deltaRadius] = holeCheck(workingBlob, curROI.Center, info.minRadius, maskCurrent);
                cents(counter,1:2) = curROI.Center+deltaCP;
                rads(counter, 1) = info.minRadius+deltaRadius;
                newCirMask= createMask(images.roi.Circle('Center',  cents(counter,1:2) , 'Radius', rads(counter, 1)), maskCurrent);

                workingBlob = bwareafilt(~newCirMask & workingBlob, 1);
                masks{counter} = newCirMask;

                % get rid of all of the points from the solution space that are within the new ROI
                tf = inROI(curROI,solutionSpaceC(:,1),solutionSpaceC(:,2));
                solutionSpaceC(tf, :) = [];
                solutionSpaceR(tf) = [];

                delete curROI % this is probably unneceessary?

                %adjust ROIindex to meet the exit condition of the while loop
                ROIindex = Inf;
                drawnow;
            else % ...if you haven't arrived at a good solution and havent exhausted search options, try the next ROI radius

                ROIindex = ROIindex+1;

            end

        end %</end while loop>

    elseif topORbot == 2
        %%
        pause(0.1);
        sz_c = size(solutionSpaceC, 1);
        if nextPtIdx<= sz_c
            C = solutionSpaceC(nextPtIdx, 1:2);
            R = solutionSpaceR(nextPtIdx);
        else
            display('MAYDAY!! ')
            disp(strcat('index is',num2str(nextPtIdx), 'but the max size of solutionSpace is', num2str(sz_c)));
        end
        [deltaCP, deltaRadius] = holeCheck(workingBlob, C, R , maskCurrent);

        hCirc = drawcircle('Center', C+deltaCP, 'Radius', R+deltaRadius);
        hROI = createMask(hCirc, maskCurrent);
        cents(counter,1:2) = hCirc.Center;
        rads(counter, 1) = hCirc.Radius;

        workingBlob = bwareafilt(~hROI & workingBlob, 1);
        masks{counter} = hROI;

        tf = isinterior(polybuffer(hCirc.Center, 'points', hCirc.Radius),solutionSpaceC(:,1),solutionSpaceC(:,2));
        solutionSpaceC(tf, :) = [];
        solutionSpaceR(tf) = [];

        delete hCirc

        drawnow;

    end
    if isempty(solutionSpaceR) || isempty(solutionSpaceC)
        return
    end
    if counter>1
        if round(cents(counter,1:2),3) == round(cents(counter-1,1:2),3)
            cents(end,:) = [];
            rads(end,:) = [];
            masks(counter) = [];
            return
        end
    end

    counter = counter+1;

    close all force
    drawnow;
end
end %/function end
%==========================================================================
% ==========================================================================


% ==========================================================================
% ==========================================================================

function [deltaCP, deltaRadius] = holeCheck(oldBlob, cp, rd, binaryMask)

newCircle = createMask(images.roi.Circle('Radius', rd, 'Center', cp), binaryMask);
newBlob = bwareafilt(~newCircle & oldBlob, 1);

filledBlob = imfill(newBlob, 'holes');
% if the filled in image is identical to the unfilled image, then there were never any holes!
containsHoles_TF = ~isequal(filledBlob, newBlob);
if ~containsHoles_TF
    deltaRadius = 0;
    deltaCP = [0,0];
    return;
end
rookMoves = [1, 2; 2, 1; 1, -2; 2, -1;...
    -1, 2; -2, 1; -1, -2; -2, -1];
noCP = [0,0];
radiusMoves = [0.5,1,4,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4];
cpMoves = [noCP; noCP; noCP;rookMoves;rookMoves;rookMoves;rookMoves];

for moveNum = 1:35
    newerC = createMask(images.roi.Circle('Radius', rd+radiusMoves(moveNum), 'Center', cp+cpMoves(moveNum, :)), binaryMask);
    newerBlob = bwareafilt(~newerC & oldBlob, 1);
    filledBlob = imfill(newerBlob, 'holes');
    containsHoles_TF = ~isequal(filledBlob, newerBlob);
    if ~containsHoles_TF
        deltaRadius = radiusMoves(moveNum);
        deltaCP = cpMoves(moveNum, :);
        disp(num2str(moveNum));
        return
    end
end
disp('unable to wiggle our way out of this hole, sire!!!')
end

% ==========================================================================
% ==========================================================================

function [nodeID, alter1, alter2] = select3Circles2reposition(circlegraph, currentCenters)

smallNetworks = find(degree(circlegraph)==1 |  degree(circlegraph)==2);
nCircles = numnodes(circlegraph)-1;

noOutBoundG = rmnode(circlegraph, nCircles+1);

% Cirlces touching the border+2cirlces are more valuable than circles with 3 alters, which we
% make clear here.
yesBounds_threeAlters = find(degree(circlegraph)==3);
noOutBounds_threeAlters = find(degree(noOutBoundG)==3);
goodThrees = intersect(yesBounds_threeAlters,noOutBounds_threeAlters );

if numel(smallNetworks)<2
    smallNetworks = [smallNetworks; goodThrees];
end

% select at random one of the these nodes with few alters.
smallNet_idxList = randperm(numel(smallNetworks));

nodeID = smallNetworks(smallNet_idxList(1));

altersOfTarget = neighbors(circlegraph, nodeID);
nAlts = numel(altersOfTarget)-1; % -1 gives you the number of alters that Aren't nodeID

switch nAlts

    case 0 %i.e. alter 2 will not be touching alter1, which at most is touching nodeID and maybe boundary
        centerP = currentCenters(nodeID, :);
        % we have already selected nodeID and Alter1. but alter1 has no alters. So we find that
        % 3rd node by minimizing the sum distance to both nodeID and alter 1
        [dxy, jj] = pdist2(currentCenters, centerP, 'euclidean', 'Smallest', 4);
        [~, minId] = min([dxy(2)+dxy(3), dxy(2)+dxy(4), dxy(3)+dxy(4)]);
        if minId <3
            alter1 = jj(2);
        else
            alter1 = jj(3);
        end
        if minId>1
            alter2 = jj(4);
        else
            alter2 = jj(3);
        end

    case 1 %i.e. we selected alter 1 randomly already. This alter1 has only 1 connect, which may or may not be outer boundary
        altersOfTarget(altersOfTarget>nCircles | altersOfTarget==nodeID) = [];
        alt1 = altersOfTarget(1);

        candidates4alter2 = neighbors(circlegraph, alt1);
        candidates4alter2(candidates4alter2==nodeID | candidates4alter2==(1+nCircles)) = [];

        if size(candidates4alter2, 1) == 1
            [x0, y0] = deal(currentCenters(alt1, 1), currentCenters(alt1, 2));
            [x1, y1] = deal(currentCenters(nodeID, 1), currentCenters(nodeID, 2));
            [x2, y2] = deal(currentCenters(candidates4alter2, 1), currentCenters(candidates4alter2, 2));
            angle_RAD = atan2(y1 - y0, x1 - x0) - atan2(y2 - y0, x2 - x0);
            angle = rad2deg(angle_RAD);
            if angle >155 && angle<205 && numel(smallNetworks)>1
                nodeID = smallNetworks(smallNet_idxList(2));
                centerP = currentCenters(nodeID, :);
                [dxy, jj] = pdist2(currentCenters, centerP, 'euclidean', 'Smallest', 4);
                [~, minId] = min([dxy(2)+dxy(3), dxy(2)+dxy(4), dxy(3)+dxy(4)]);
                alter1 = jj(2);
                alter2 = jj(3); %(jj(1)=nodeID)
                if minId==2
                    alter2 = jj(4);
                elseif minId == 3
                    alter1 =jj(4);
                end
            else
                alter1 = alt1;
                alter2 = candidates4alter2(1);
            end

        elseif size(candidates4alter2, 1) == 0 % alter1 touches only the boundary and nodeID only
            CPs = [currentCenters(candidates4alter2(:),1:2) ; currentCenters(nodeID, 1:2)];
            [~, gg] = pdist2(CPs, currentCenters, 'euclidean', 'Smallest', 1);
            alter1 = setdiff(gg(2:3),[nodeID, altersOfTarget]);
            alter2 = altersOfTarget(1);

        else
            [~, idx] = pdist2(currentCenters(candidates4alter2, :), currentCenters(nodeID, :), 'euclidean', 'Smallest', 2);
            alter1 = altersOfTarget(1);
            alter2 = candidates4alter2(idx(1));
        end

    case 2
        altersOfTarget(altersOfTarget>nCircles | altersOfTarget==nodeID)=[];
        alter1 = altersOfTarget(1);
        alter2 = altersOfTarget(2);

    otherwise % >2
        disp('not there yet');

%
%
%
%
%

end

assert(nodeID~=alter1);assert(nodeID~=alter2);
end

% ==========================================================================
% ==========================================================================

function [triangulationRing, orderIDs_cw] = makeLocatingRing(centerPoint, radius)
% makes a donut shape with 48 triangles 
orderIDs_cw = zeros(3, 48,'double');

% precompute some things that are invariant throughout each loop iteration
    p = linspace(0, 2.*pi,25);
    x = (radius+400)*cos(p)'+centerPoint(1,1);
    y = (radius+400)*sin(p)'+centerPoint(1, 2);
    x2 = (radius/2)*cos(p)'+centerPoint(1,1);
    y2 = (radius/2)*sin(p)'+centerPoint(1,2);

for k = 1:3
   rotDeg = (k-1)*-5;
    ps = rotate(polyshape(x,y), 22.5+rotDeg, centerPoint(1, 1:2));
    ps2 = rotate(polyshape(x2,y2), rotDeg, centerPoint(1, 1:2));
    
    a = (1:1:24)';b  = (25:1:48)';
    constraints = [a, circshift(a, -1,1); b, circshift(b, -1,1)];
    
    DT = delaunayTriangulation([ps.Vertices; ps2.Vertices], constraints);
    inside = isInterior(DT);
    tr = triangulation(DT(inside,:),DT.Points);
    triangulationRing{k} = tr;
    
    %% now we can recursively ascertain the order of the IDs for each triang. 

    % Recursion requires requires 2 known values on which to base next value (the starting
    % point and then the direction CW or CCW). maybe this would be more readable if I did n-2,
    % n-1, and n.... but I did n-1, n, n+1. sorry
    
    orderIDs_cw(k, 1) = pointLocation(tr, centerPoint(1, 1)-radius-350, centerPoint(1, 2)+60);
    % notice we didn't chose a central point here. Rather, this point is located squarely in
    % the intersection of what should be the first triangle rotated 0,-5,-10 deegrees (i.e.
    % CW direction), such that we can meaningfully combine these 3 arrays later. 
    N = neighbors(tr, orderIDs_cw(k, 1));
    N = N(~isnan(N));
    incents = incenter(tr, N');
    
    % So I take a point due west of the center point (just farther than the inner radius so its
    % "IN" sans doubte. then I find its trigangle. this will be new #1. then, I have 2 IDs to
    % choose from, 1 up and 1 down. So I query their y coords and pick the bigger one to be 2nd.
    [~, maxId] = max([incents(1, 2), incents(2,2)]);
    orderIDs_cw(k, 2) = N(maxId);
    
    % this is my recursive loop to generate the other 46 triangle IDs
    for counter = 2:47
    N = neighbors(tr, orderIDs_cw(k, counter));
    N = N(~isnan(N)); % bc of the constrains one edge will always be NaN
    orderIDs_cw(k, counter+1) = N(N~=orderIDs_cw(k, counter-1)); 
    % this says: don't repeat the last triangle, ie if it cant be the previous or the nan, so by 
    % process of elimination we get our next triangle 
    end
    
end




% the angles the rays make from center to the centroid is
% %angles = (187.5:-7.5:-165)';
end


% ==========================================================================
% ==========================================================================

function direction = queryBlobs(masks, centers, radii, binMask_blob, binMask_tissue, info )
% direction = query3Blobs(centers, radii, info, binaryMask)

% loop through each of the new circles you've just added to the packing
    for k= 1:size(radii,1)

        cp = centers(k, 1:2);
        radius = radii(k,1);

        % makeLocationRing is a triangulation we will use to hone in on our target direction
        [allTriangulations, orderedIDs_cw] = makeLocatingRing(cp, radius);

        % define a region of the image larger than 
       searchArea = createMask(images.roi.Circle('Radius', radius+410, 'Center', cp), binMask_tissue);
       [rws, cols] = find(binMask_blob & searchArea);
       
       combinedData = zeros(72, 1, 'double');

       % loop through each of the 3 slightly offset triangulations in order to combine them
       for tri = 1:3

       ID = pointLocation(allTriangulations{tri}, cols, rws);
    
       % identify which simplex is the most popular "closest" point of the triuangulation using histcounts. 
       [n, edges] = histcounts(ID, 48);
        outerIDs = orderedIDs_cw(tri, 1:2:end);
       orderedBinCounts = n(outerIDs);
       data_w_zeros = reshape([orderedBinCounts; zeros(2, 24)], [], 1);
        if tri == 1
        combinedData = data_w_zeros;
        else
       combinedData = combinedData+circshift(data_w_zeros, tri-1);
        end
       end

       %% NOTE ALL OF THESE X VALUES ARE OFF BY -2.5, NEED TO ADD THAT BEFORE USING
% Create the polar histogram
signal = flipud(combinedData);

% pad signal - I did it this a silly way just because its late...
padSignal = [signal; signal; signal];
xVals = (0:5:719)'; % this is the centers, not the edges, which need a 2.5deg offset

% clip pad Signal to match xVals
padSignal = padSignal(36: 179);

% convolve filter and de-pad
F = griddedInterpolant(xVals,padSignal,  'pchip');
xq = linspace(0,718, 719);
interpValues = F(xq);
interp_noPad = interpValues(:, 181:540);
interp_noPad_xCenters = (180:1:539)';

[indMin,PromMn] = islocalmin(interp_noPad);
minimaLocations = interp_noPad_xCenters(find(indMin)); 
minimaLocations0_360  = mod(minimaLocations, 360);
minimaValues = interp_noPad(indMin);
MINprominence = PromMn(indMin);
kissPointsAlters = minimaLocations(MINprominence>10000);
potentialAlters2Kiss = minimaLocations(MINprominence<=10000);

[indMx, PromMx] = islocalmax(interp_noPad);
maxes = mod(interp_noPad_xCenters(find(indMx)), 360);

%if numel(potentialAlters2Kiss) == 3 || numel(potentialAlters2Kiss) == 4
queryMinima = mod(sort(potentialAlters2Kiss,'ascend'), 360);

% what I'm doing here is trying to identify one of the easiest repositioning targets, two
% adjacent cirlces I can reposition circle to touch both of
for p = 1:(numel(potentialAlters2Kiss))
    if p < numel(potentialAlters2Kiss)
        maxInRange = find(maxes>queryMinima(p,1) & maxes<queryMinima(p+1,1));
    else
        maxInRange = [find(maxes>queryMinima(p,1) ), find(maxes<queryMinima(1,1))];
    end
    nMax(p,1) = numel(maxInRange);
    maxesInRange_cell{p} = maxInRange;
end
id1max = find(nMax(:) == 1);
id1 = maxesInRange_cell{id1max};
centerpeak = maxes(id1);
direction = [queryMinima(id1max); centerpeak;queryMinima(id1max+1)];

Update(binMask_blob, triangulation1, interp_noPad, direction)

%% the next steps are to:
% move in the directio of Direction vector enough to touch both those cirles. 
% then, expand radius to at least touch a third 
    
       

    


    end
end

% ==========================================================================
% ==========================================================================

function plotUpdate(binMask_blob, triangulation1, interp_noPad, direction)
triangulation1 = allTriangulations{1};
%plot it as a polar histogram
[x,y] = pol2cart(deg2rad(direction+2.5), [250;500;250]);
centerTri = sum(triangulation1.Points, 1)./size(triangulation1.Points, 1);
x = x+centerTri(:,1); y = y+centerTri(:,2);
psDir = polyshape(x,y);

f = figure;
ax1 = axes(f);
imshow(binMask_blob, 'Parent', ax1); hold on; 
triplot(allTriangulations{1}, 'Parent', ax1);
plot(psDir);


hold off
subplot(1,2,1, ax1)

pax =subplot(1,2,2, polaraxes);
polarhistogram(pax,'BinEdges', pi+deg2rad(-0.5:1:359.5), 'BinCounts', interp_noPad);
pax.ThetaDir = 'clockwise'; drawnow;

end

% ==========================================================================
% ==========================================================================




