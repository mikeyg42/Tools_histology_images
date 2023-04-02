function [s, derivM, derivF] = curveFittingOfTissueBorders(pointsPerSide,cornerN, PointsToFit_move, PointsToFit_fix, movingIm, fixedIm)
%syntax: s = curveFittingOfTissueBorders(pointsPerSide, cornerN,PointsToFit_move, PointsToFit_fix, movingIm, fixedIm) 

if mod(cornerN, 2)==0
    % this step sets the y values to be the x value for the left and right sides of the tissue, in effect rotating the image 90ª
    PointsToFit_move = fliplr(PointsToFit_move);
    PointsToFit_fix= fliplr(PointsToFit_fix);
end

%% we've arrived at the interpolation station
% first step we rotate the edge points so that our two corners are zero crossings .... this will make
% our lives easier later on

% ================== first, rotate MOVING point cloud ==================
% 1. define 2 vectors to find angle between and define subsequent theta of rotation with. 
PointsShifted_move = PointsToFit_move(:, 1:2) - [0, PointsToFit_move(end, 2)];
toBeMovedToZero_m = PointsShifted_move(1, 1:2);
pivotM = PointsShifted_move(end, 1:2);

% 2. calculate theta
vec1 = pivotM - toBeMovedToZero_m;
vec2 = [toBeMovedToZero_m(1), 0] - pivotM;
theta = atan2d(norm(cross([vec1, 0],[vec2, 0])),dot(vec1,vec2)); %gives the CLOCKWISE rotation

% 3. rotate first moving the axis of rotation to "pivot" point. After this transformation, pivot should
%   be the only point in the point cloud unmoved and first and last points should have
%   equal y
PointsRotated_move = rotate_CCW_around_a_point(PointsToFit_move, -theta, PointsToFit_move(end, 1:2));

% -------... brief discussion about Algebra, OPTIONAL read ....<start discussion> ----------
% topic: why bother with Chebychev nodes? 

% Of Pafnuty Chebyshev's prolific contributions to math (e.g.
% top five most famous inequalities "Chebyshev inequality", his work with prime numbers,  probability...)
% were the so-called Chebyshev polynomials.
% They all exist within [-1, 1], and are orthogonal/
% For "the first kind" of this class of polynomials, the x-coordinates of the roots of the
% n'th polynomial in this series (where T sub-n of (cos(theta)) = cos(n*theta)) elegantly 
% work out to be n equally spaced points along a unit semi-circle . 
% 
% Just imagine placing evenly subdividing the rim of half a circle and 
% then looking at those points' x-coordinates. Near -1 and 1 the points will be tightly
% packed, but on the curve they will be equidistant. 
%       Chebyshev points:
%        _*--^^*^^--*_   
%      */ |    |    | \*  <--- the asterix in theory equally spaced apart (use imagination)
%     *|  |    |    |  |*        
%     oo--o----o----o--oo <--- "diameter" of circle with projected points irregularly spaced apart (nodes = no's)
%     
% Why do this? There is a phenomenon described by a german mathematician named 
% Carl David Tolme Runge in 1901 (Runge's phenomenon) where one's efforts at
% interpolating a high degree polynomial with points equidistant along the x-axis are 
% thwarted on the edges by large, expanding oscillations. Chebyshev nodes end up being one way to guarantee
% this does NOT happen. As a result, they're commonly used for interpolating high degree
% polynomials. 

% to calculate n number of Chebyshev nodes, each point k in range [1:1:n] is:
% x(k) = cos((2k-1)*(pi/2n). I use the nice matlab function cospi((2k-1)/2n)= x for
% each n! 
% ---------------------- ** </end discussion> ** -----------------------


xcoords_movingCheby = ChebyshevNodes(PointsRotated_move(1,1), PointsRotated_move(end, 1), pointsPerSide);
%note, ChebyshevNodes will generate n-1 nodes, and then will add in the corner points,
%results in n+1 nodes. Later one, we will remove 1 of the 2 corners from each line so that
%we have one corner per side and nothing is double counted

% 4. fit rotated points to a curve, trying:
%     { 1->3   1st-3rd polynomial of the 1st, 2nd, & 3rd degrees, 
%       4      4th degree polynomial with fixed roots at min/max of x,
%       5      sine wave. }

%next we fit the curve to the entirety of edges data (1st step of FITTHECURVE is a subsampling...)
[curveFcn_M, derivM] = fitTheCurve(PointsRotated_move(:, 1), PointsRotated_move(:, 2), cornerN, "Moving");

% 5. use curve to project the Chebychev nodes and find estimated y
try y_Cheb_m = polyval(curveFcn_M, xcoords_movingCheby(:)); %for polynomial curves
catch
   y_Cheb_m = feval(curveFcn_M, xcoords_movingCheby(:)); %for nonpolynomials
end
% 6. rotate points around same pivot back to where they belong (CCW), before any rotation
xyPoints_M = rotate_CCW_around_a_point([xcoords_movingCheby', y_Cheb_m], theta, PointsToFit_move(end, 1:2));

% ================== 2nd, repeat with fixed ==================
PointsShifted_fix = PointsToFit_fix(:, 1:2) - [0, PointsToFit_fix(end, 2)];
toBeMovedToZero_f = PointsShifted_fix(1, 1:2);
pivotF = PointsShifted_fix(end, 1:2);

vec3 = pivotF - toBeMovedToZero_f;
vec4 = [toBeMovedToZero_f(1), 0] - pivotF;
theta2 = atan2d(norm(cross([vec3, 0],[vec4, 0])),dot(vec3,vec4)); %gives the CLOCKWISE rotation

PointsRotated_fix = rotate_CCW_around_a_point(PointsToFit_fix, -theta2, PointsToFit_fix(end, 1:2));

xcoords_fixedCheby = ChebyshevNodes(PointsRotated_fix(1,1), PointsRotated_fix(end, 1), pointsPerSide);

[curveFcn_F, derivF] = fitTheCurve(PointsRotated_fix(:, 1), PointsRotated_fix(:, 2), cornerN, "Fixed");

try y_Cheb_f = polyval(curveFcn_F, xcoords_fixedCheby(:)); %for polynomials
catch
y_Cheb_f = feval(curveFcn_F, xcoords_fixedCheby(:)); %for nonpolynomials
end

xyPoints_F = rotate_CCW_around_a_point([xcoords_fixedCheby', y_Cheb_f], theta2, PointsToFit_fix(end, 1:2));

%===================================================================%


%============= Force Points Back on to bounary of masks ============%
    % PROBLEM: despite the relative orientation of our control points being well set,
    % they've inevitable DRIFTED away from the edge of our tissue mask.
    % SOLUTION: force each edge point back to the mask edge, while preserving 
    % their relative orietnation.
    % PROCEDURE: we can use relative arc lengths
    
    FarcLengths = zeros(size(xyPoints_F,1),1);
    MarcLengths = zeros(size(xyPoints_M,1),1);
    
    % I had to put this into a try-catch block because the polynomial functions work much
    % better in this form but the sine curve I need to keep as an anonymous function. it
    % "catches" only when there the fit curve isn't a polynomial
    
    %for p = 2:size(xyPoints_F,1)
    %  MarcLengths(p,1) = arcLength(derivM, xyPoints_M(1,1), xyPoints_M(2,1));
    % MarcLengths(p,1) = arcLength(derivM, min(xyPoints_M(1,1), xyPoints_M(p,1)), max(xyPoints_M(1,1), xyPoints_M(p,1)));
    % FarcLengths(p,1) = arcLength(derivF, xyPoints_F(p,1), xyPoints_F(p-1,1));
    %end
    fdiff = xyPoints_F(:,1:2) - circshift(xyPoints_F(:,1:2),1,1);
    mdiff = xyPoints_M(:,1:2) - circshift(xyPoints_M(:,1:2),1,1);
    
    Farcsegments = hypot(fdiff(2:end,1), fdiff(2:end,2));
    Marcsegments = hypot(mdiff(2:end,1), mdiff(2:end,2));
    
    fullArcLength_F = sum(Farcsegments(:));
    fullArcLength_M = sum(Marcsegments(:));
    
    deltaIdx = find(Marcsegments < 7 | Farcsegments < 7); %arbitrarily chosen theshhold
    
    %we can't lose our corner points. so here we index the other point contributing to the tiny distance
    deltaIdx(deltaIdx == 1) = 2; 
    deltaIdx(deltaIdx == size(Marcsegments,1)) = size(Marcsegments,1)-1;
    
    %convert from segments to lengths with cumsum
    FarcLengths = cumsum(Farcsegments);
    MarcLengths = cumsum(Marcsegments);
    
    % add a zero to the top of this list we lost with the circshift line. We could not add
    % til now because the 0 distance would have been under threshold
    FarcLengths = [0; FarcLengths]; 
    MarcLengths = [0; MarcLengths]; 
    deltaIdx = deltaIdx+1; % accomodating the added 0 by shifting all index up
    
    fratios = FarcLengths./fullArcLength_F;
    mratios = MarcLengths./fullArcLength_M;
    
    % Use the ratios of segment lengths to transpose from curves into the real data 
    fidx = round(fratios.*size(PointsToFit_fix,1)); fidx(1,1) = 1;
    midx = round(mratios.*size(PointsToFit_move,1)); midx(1,1) = 1;
    if midx(2)<5  || fidx(2)< 5
        deltaIdx = [2, deltaIdx];
    end
    
    %we take the upper indicies that are too close to neighbors, and wiggle down by taking
    %2/3rd 1/3rd weighted average with the next index
    topDelta = deltaIdx>size(Marcsegments,1)/2;
    fidx(deltaIdx(topDelta)) = ceil((2*fidx(deltaIdx(topDelta))+fidx(deltaIdx(topDelta)-1))/3);
    midx(deltaIdx(topDelta)) = ceil((2*midx(deltaIdx(topDelta))+midx(deltaIdx(topDelta)-1))/3);
    
    %we take the lower indices and move them up closer to the next point (by 1/3)
    lowerDelta = deltaIdx<size(Marcsegments,1)/2;
    fidx(deltaIdx(lowerDelta)) = floor(( 2*fidx(deltaIdx(lowerDelta))+fidx(deltaIdx(lowerDelta)+1) )/3 );
    midx(deltaIdx(lowerDelta)) = floor(( 2*midx(deltaIdx(lowerDelta))+midx(deltaIdx(lowerDelta)+1) )/3 );
    
    newPointsF = PointsToFit_fix(fidx,:); %these are all XY not RC
    newPointsM = PointsToFit_move(midx,:);
    
    newPointsM = cpcorr(newPointsM, newPointsF, movingIm, fixedIm);
    newPointsF = cpcorr(newPointsF, newPointsM, fixedIm, movingIm);
    
%============ now we determine middle grid points ==================%
    if mod(pointsPerSide, 2)==0 %then nPoints+1 is odd:
        fMid = 1+(pointsPerSide/2);
        mMid = 1+(pointsPerSide/2);
    else
        opts = [(pointsPerSide+1)/2,(pointsPerSide+1)/2+1]; % nPoints+1 is even
        pointoptions = [newPointsF(opts, :) - newPointsM(opts, :);newPointsF(opts, :) - flipud(newPointsM(opts, :))];
        [~, ix] = min(hypot(pointoptions(:, 1), pointoptions(:, 2)));
        if ix == 1 || ix == 3, fMid = opts(1); else, fMid = opts(2); end
        if ix == 1 || ix == 4, mMid = opts(1); else, fMid = opts(2); end
    end
    m_Ind(1) = floor((1+mMid)/2); f_Ind(1) = floor((1+fMid)/2);
    m_Ind(2) = ceil((pointsPerSide+1+m_Ind)/2);f_Ind(2) = ceil((pointsPerSide+1+f_Ind)/2);
    
    midGridf(1, 1:2) = newPointsF(f_Ind(1), 1:2);
    midGridf(2, 1:2) = newPointsF(fMid, 1:2);
    midGridf(3, 1:2) = newPointsF(f_Ind(2), 1:2);
    
    midGridm(1, 1:2) = newPointsM(m_Ind(1), 1:2);
    midGridm(2, 1:2) = newPointsM(fMid, 1:2);
    midGridm(3, 1:2) = newPointsM(m_Ind(2), 1:2);

%================ DONE w/ setting Points!! =======================%


%=========final step: FLIP XY BACK (on 2nd and 4th sides) ==========%
% need to flip the x and the y coordinates back for the two sides 
% where we earlier swapped x and  (the second and fourth)

if mod(cornerN, 2)==0 
    newPointsM = fliplr(newPointsM);
    newPointsF = fliplr(newPointsF);
    midGridm = fliplr(midGridm);
    midGridf = fliplr(midGridf);
end
%===================================================================%


%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%
%% saving the results into a structural array "s"

s = struct('ID', [], 'middleGrid', [], 'xyPoints', []);

s(1).ID = 'mixed';
s(1).middleGrid = midGridm;
s(1).xyPoints = newPointsM;

s(2).ID = 'fixed';
s(2).middleGrid = midGridf;
s(2).xyPoints = newPointsF;

%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%

end


function middleGridPoints = setMiddleGridPoints(final_xyPoints, deriv_fcn) %(x,y)
%Arc Length!!
% calculating dy/dx numerically and integrating using trapz
fullArcLen = arcLength(deriv_fcn, final_xyPoints(1,1),final_xyPoints(end,1));
    
sz = size(final_xyPoints,1);
halfIdx = floor(sz/2);


arcLen1 = arcLength(deriv_fcn, final_xyPoints(1,1),final_xyPoints(halfIdx,1));
arcLen2 = arcLength(deriv_fcn, final_xyPoints(halfIdx,1),final_xyPoints(end,1));

if arcLen1 < fullArcLen/2
   flg=1;
else
   flg=2;
end

k=0; %basically we keep oscillating here until we have setteled on the place to split the arc's that is most even
while (flg==1 && arcLen1 < fullArcLen/2 && arcLen2 > fullArcLen/2) || (flg==2 && arcLen1 > fullArcLen/2 && arcLen2 < fullArcLen/2)
    if flg==1
        k = k+1;
    else
        k = k-1;
    end
    idx = halfIdx+k*flg;
    arcLen1 = arcLength(deriv_fcn, final_xyPoints(1,1),final_xyPoints(idx,1));
    arcLen2 = arcLength(deriv_fcn, final_xyPoints(idx,1),final_xyPoints(end,1));
end

%note these points, being averages of the midpoint and end points on an arc, will be OFF
%the arc considerably. That is why our next function call is to the pdist2 to locate the
%closest point ON the arc to situate our quarter points!
quarterPoints = [round((final_xyPoints(1,:)+final_xyPoints(halfIdx+k*flg, :))./2);... 
                 round((final_xyPoints(halfIdx+k*flg,:)+final_xyPoints(end, :))./2)];
             
[~, ind] = pdist2(final_xyPoints,quarterPoints,'euclidean', 'Smallest', 1);
             
middleGridPoints = zeros(3, 2, 'double');
middleGridPoints(1, 1:2) = final_xyPoints(ind(1), 1:2);
middleGridPoints(2, 1:2) = final_xyPoints(halfIdx+k*flg, 1:2);
middleGridPoints(3, 1:2) = final_xyPoints(ind(2), 1:2);

arcLen4 = arcLength(deriv_fcn, middleGridPoints(1, 1),middleGridPoints(2, 1));
arcLen5 = arcLength(deriv_fcn, middleGridPoints(2, 1),middleGridPoints(3, 1));
if max(arcLen4, arcLen5)/(arcLen4 + arcLen5) > 0.7 || (arcLen4 + arcLen5)/fullArcLen >0.7 || arcLen4/arcLen1 >0.7 || arcLen5/arcLen2 > 0.7
    %^ all of these ratios SHOULD be 0.5
    disp('middle arc Lengths are very much NOT balanced evenly ....');
    pause(10);
end
end


