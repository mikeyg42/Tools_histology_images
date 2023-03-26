function s = curveFittingOfTissueBorders(pointsPerSide,cornerN, PointsToFit_move, PointsToFit_fix, movingIm, fixedIm)
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

% -------...intermission: brief lecture on Real Algebra....class is now IN SESSION ----------
% instead of= using linspace to sample, its been an helpful to use Chebyshev points for nodes!
% Imagine drawing a circle s.t. the two points are the diameter. 
% Of Chebychev's many brilliant contributions to algebra was the Chebychev polynomials.
% For "the first kind" of this class of polynomials, the x-coordinates of the roots of the
% n'th polynomial in this series (where T sub-n of (cos(theta)) = cos(n*theta)) elegantly 
% work out to be n equally spaced points along a unit semi-circle . 
% Confused? Sorry! Just imagine placing evenly subdividing the rim of half a circle and 
% then looking at those points' x-coordinates. Near -1 and 1 the points will be tightly
% packed, but on the curve they will be equidistant. 
%       Chebychev points:
%        _*--^^*^^--*_   
%      */ |    |    | \*  <--- the asterix in theory equally spaced apart (use imagination)
%     *|  |    |    |  |*        
%     oo--o----o----o--oo <--- "diameter" of circle with projected points irregularly spaced apart (nodes = no's)
%     
% Why do this? There is a phenomenon described by a german mathematician named 
% Carl David Tolme Runge in 1901 (Runge's phenomenon) where one's efforts at
% interpolating a high degree polynomial with points equidistant along the x-axis are 
% thwarted on the edges by huge, expanding oscillations. Chebychev nodes oddly enough are one way to guarantee
% this doesnt happen. As a result, they're commonly used for interpolating (which is what
% we are doing) high degree polynomials (a reasonable way to describe our tissue borders.
% to calculate n number of points, each point k in range [1:1:n] is:
% x(k) = cos((2k-1)*(pi/2n). I use the nice matlab function cospi((2k-1)/2n)= x for
% each n! 
% ---------------------- **bell rings, class is over** ----------------------------

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
y_Cheb_m = feval(curveFcn_M, xcoords_movingCheby);

% 6. rotate points around same pivot back to where they belong (CCW), before any rotation
final_Points_M = rotate_CCW_around_a_point([xcoords_movingCheby; y_Cheb_m]', theta, PointsToFit_move(end, 1:2));

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
y_Cheb_f = feval(curveFcn_F, xcoords_fixedCheby);
final_Points_F = rotate_CCW_around_a_point([xcoords_fixedCheby; y_Cheb_f]', theta2, PointsToFit_fix(end, 1:2));

%============================================================%
% need to flip the x and the y coordinates back for the two sides 
% where we earlier swapped x and y
xyPoints_m = final_Points_M;
xyPoints_f = final_Points_F;

if mod(cornerN, 2)==0 
    xyPoints_m = fliplr(xyPoints_m);
    xyPoints_f = fliplr(xyPoints_f);
end

%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%
%% saving the results into a structural array "s"

s = struct('ID', [], 'middleGrid', [], 'xyPoints', []);

s(1).ID = 'middle';
s(1).middleGrid = setMiddleGridPoints(final_Points_M, derivM);
s(1).xyPoints = xyPoints_m;

s(2).ID = 'fixed';
s(2).middleGrid = setMiddleGridPoints(final_Points_F, derivF); 
s(2).xyPoints = xyPoints_f;

%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%

end


function middleGridPoints = setMiddleGridPoints(final_Points, deriv) %(x,y)
%Arc Length!!
% calculating dy/dx numerically and integrating using trapz
fullArcLen = arcLength(deriv, final_Points(1,1),final_Points(end,1));
    
sz = size(final_Points,1);
halfIdx = floor(sz/2);

arcLen1 = arcLength(deriv, final_Points(1:halfIdx,1),final_Points(1:halfIdx,2));
arcLen2 = arcLength(deriv, final_Points(halfIdx:end,1),final_Points(halfIdx:end,2));

if arcLen1 < fullArcLen/2
   flg=1;
else
   flg=2;
end

k=0;
while (flg==1 && arcLen1 < fullArcLen/2 && arcLen2 > fullArcLen/2) || (flg==2 && arcLen1 > fullArcLen/2 && arcLen2 < fullArcLen/2)
    if flg==1
        k = k+1;
    else
        k = k-1;
    end
    arcLen1 = arcLength(deriv, final_Points(1:halfIdx+k*flg,1),final_Points(1:halfIdx+k*flg,2));
    arcLen2 = arcLength(deriv, final_Points(halfIdx+k*flg:sz,1),final_Points(halfIdx+k*flg:sz,2));
end


quarterMark = round(Index_crossing/2); 
threequartersMark = round((Index_crossing + sz)/2);

middleGridPoints = zeros(3, 2, 'double');
middleGridPoints(1, 1:2) = final_Points(quarterMark, 1:2);
middleGridPoints(2, 1:2) = final_Points(Index_crossing, 1:2);
middleGridPoints(3, 1:2) = final_Points(threequartersMark, 1:2);
end