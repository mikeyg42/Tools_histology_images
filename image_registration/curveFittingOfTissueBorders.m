function s = curveFittingOfTissueBorders(pointsPerSide, cornerNum, movingM, fixedM, PointsToFit, PointsToFit_fix, varargin)
%syntax: s = curveFittingOfTissueBorders(pointsPerSide, cornerNum, movingM, fixedM,...
%                   PointsToFit, PointsToFit_fix, varargin)
% OR:    s = curveFittingOfTissueBorders(pointsPerSide, cornerNum, movingM, fixedM, ...
%                   PointsToFit, PointsToFit_fix, IMG_gray, LESION_gray)

%varargin is mainly for visualization purposes

try
    if ismatrix(varargin{1}) && ismatrix(varargin{2})
        IMG_gray = varargin{1};
        LESION_gray = varargin{2};
    end
catch
end

cornerN_Plus1 = cornerNum + 1;
if cornerN_Plus1==5
    cornerN_Plus1 = 1;
end

mStart = movingM(cornerNum, 1:2);
mFinish = movingM(cornerN_Plus1, 1:2);
fStart = fixedM(cornerNum, 1:2);
fFinish = fixedM(cornerN_Plus1, 1:2);

if mod(cornerNum, 2)==0
    % this step sets the y values to be the x value for the left and right sides of the tissue, in effect rotating the image 90ª
    mStart = fliplr(mStart);
    mFinish = fliplr(mFinish);
    fStart = fliplr(fStart);
    fFinish = fliplr(fFinish);
    PointsToFit = fliplr(PointsToFit);
    PointsToFit_fix= fliplr(PointsToFit_fix);
end

% instead of= using linspace, I wrote a fcn to extract the Chebyshev
% nodes between two points. These nodes are result of drawing a circle s.t. the two points are the diameter. 
% Mr.Chebychev placed some points along the rim of the circle and 
% then projected them down to the diameter. 
%       Chebychev points:
%        _*--^^*^^--*_   
%      */ |    |    | \*  <--- the asterix in theory equally spaced apart (use imagination)
%     *|  |    |    |  |*        
%     oo--o----o----o--oo <--- "diameter" of circle with projected points irregularly spaced apart (nodes = no's)
%     
% they are more concentrated near the edges and have fewer in the middle. 
% given we wil be projecting back out onto a curve very much like this, it seems a propos.

qPoints_m = ChebyshevNodes(mStart(:,1), mFinish(:, 1), pointsPerSide+1);
% if sign(mFinish-mStart)~=sign(PointsToFit(end, 1)-PointsToFit(1, 1)) %ie if X is increasing but PointsToFit aint
%     qPoints_m = flipud(qPoints_m);
% end

qPoints_f = ChebyshevNodes(fStart(:, 1), fFinish(:, 1), pointsPerSide+1);
% if sign(fFinish-fStart)~=sign(PointsToFit(end, 1)-PointsToFit(1, 1)) %ie if X is increasing but PointsToFit aint
%     qPoints_f = flipud(qPoints_f);
% end

%% curve fitting time
% first we rotate the image so that our two corners are symmetric (both can be a root) .... 
% We want the Corner Points to be part of the curve. and best way to do
% that is to insist that corners are "zero-crossings" ie roots.

% ================== first, transform moving ims ==================
% 1. define 2 vectors to find angle between and define subsequent theta of rotation with. 
pivot = PointsToFit(end, 1:2);

% 2. calculate theta
toBeSetToZero = PointsToFit(1, 1:2);
vec1 = toBeSetToZero - pivot;
vec2 = [toBeSetToZero(1), 0];
theta = acosd(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
if theta>90
    theta = 180-theta; %this will happen on the 3rd and 4th corners
end
if cornerNum==2 || cornerNum == 4
    theta = -theta;
end

% 3. rotate first moving the axis of rotation to "pivot" point. After this transformation, pivot should
%   be the only point in the point cloud unmoved
PointsToFit_norm = rotate_CCW_around_a_point(PointsToFit, theta, pivot);

% 4. fit these rotated points to the best curve {
%    1. sine wave,    2-4. polynomial of the 1st, 2nd, & 3rd degrees, 

xm = PointsToFit_norm(:, 1);
ym = PointsToFit_norm(:, 2); 

curveFcn_M = fitTheCurve(LESION_gray, xm, ym, cornerNum, "Moving");

% 5. use curve to project these points x ->y. use arc length integration
%   to ensure points are equidistant spacing.
xQueryM1 = min(xm):0.5:max(xm);
xQueryM2 = qPoints_m;

curve_Points_m = feval(curveFcn_M, xQueryM1);
Cheb_Points_m = feval(curveFcn_M, xQueryM2);

% NO CLUE WHY I NEEED THIS? something wacky is happening?
if size(Cheb_Points_m', 2) == 1
    Cheb_Points_m = [xQueryM2;Cheb_Points_m]';
end

% 6. rotate points back CW to where they belong, before you rotated them
final_Points_M = rotate_CCW_around_a_point([xQueryM1; curve_Points_m]', -theta, pivot);
Cheb_Points_m = rotate_CCW_around_a_point(Cheb_Points_m, -theta, pivot);


% ================== 2nd, repeat with fixed ==================

% 1. define 2 vectors to find angle between and define subsequent theta of rotation with. 
fpivot = PointsToFit_fix(end, 1:2);

% 2. calculate theta to rotate the domaine of points to fit
toBeZero_fixed = PointsToFit_fix(1, 1:2);
vec3 = toBeZero_fixed - fpivot;
vec4 = [toBeZero_fixed(1), 0];
theta2 = acosd(dot(vec3,vec4)/(norm(vec3)*norm(vec4)));
if theta2>90
    theta2 = 180-theta2; 
end
if cornerNum==2 || cornerNum == 4
    theta2 = -theta2;
end

% 3. rotate first moving the axis of rotation to "pivot" point. After this transformation, pivot should
% be the only point in the point cloud unmoved
PointsToFit_fixNorm = rotate_CCW_around_a_point(PointsToFit_fix, theta2, fpivot);

% 4. fit rotated points to the best curve {
%    1. sine wave,      2-4. polynomial of the 1st, 2nd, & 3rd degrees, 
xf = PointsToFit_fixNorm(:, 1);
yf = PointsToFit_fixNorm(:, 2);

curveFcn_f = fitTheCurve(IMG_gray, xf, yf, cornerNum, "Fixed");

%5. use curve to transform you chebychev points. use arc length integration
%to ensure points are equidistant spacing.
xQueryF1 = min(xf):0.5:max(xf);
xQueryF2 = qPoints_f;

curve_Points_f = feval(curveFcn_f, xQueryF1);
Cheb_Points_f = feval(curveFcn_f, xQueryF2);

if size(Cheb_Points_f', 2) == 1
    Cheb_Points_f = [xQueryF2;Cheb_Points_f]';
end

%6. rotate points back CW to where they belong, before you rotated them
final_Points_F = rotate_CCW_around_a_point([xQueryF1; curve_Points_f]', -theta2, fpivot);
Cheb_Points_f = rotate_CCW_around_a_point(Cheb_Points_f, -theta2, fpivot);


%============================================================%
% need to flip the x and the y coordinates back for the two sides 
% where we earlier swapped x and y
xyPoints_m = Cheb_Points_m;
xyPoints_f = Cheb_Points_f;

if mod(cornerNum, 2)==0 
    xyPoints_m = fliplr(xyPoints_m);
    xyPoints_f = fliplr(xyPoints_f);
end

%% VISUALIZING RESULTS

%rotated poinnts fit to the curve
fig4 = uifigure('Visible', 'on');
tl = tiledlayout(fig4, 2,2,'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
ax11 = axes(tl);
title(ax11, 'FIXED');
imshow(IMG_gray, 'Parent', ax11); 
hold on; plot(xf, yf, 'c-o', yf, xf, 'm-o'); hold off

nexttile
ax12 = axes(tl);
title(ax12, 'MOVING');
imshow(LESION_gray, 'Parent', ax12); 
hold on; plot(xm, ym, 'c-o', ym, xm, 'm-o'); hold off

nexttile
ax21 = axes(tl);
imshow(IMG_gray,'Parent', ax21); 
hold on; plot(Cheb_Points_f(:,1), Cheb_Points_f(:,2), 'go', 'MarkerSize', 10); hold off

nexttile
ax22 = axes(tl);
imshow(LESION_gray, 'Parent', ax22); 
hold on; plot(Cheb_Points_m(:,1),Cheb_Points_m(:,2), 'mo', 'MarkerSize', 10); hold off

fig4.Visible = 'on';
drawnow limitrate nocallbacks
%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%
%% saving the results into a structural array "s"

s = struct('ID', [], 'middleGrid', [], 'xyPoints', []);

s(1).ID = 'middle';
s(1).middleGrid = setMiddleGridPoints(final_Points_M);
s(1).xyPoints = xyPoints_m;

s(2).ID = 'fixed';
s(2).middleGrid = setMiddleGridPoints(final_Points_F); 
s(2).xyPoints = xyPoints_f;

%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%

end


function middleGridPoints = setMiddleGridPoints(final_Points)
%Arc Length!!
% calculating dy/dx numerically and integrating using trapz
fullArcLen = arcLength(final_Points(:,1),final_Points(:,2));
    
sz = size(final_Points,1);
halfIdx = floor(sz/2);

arcLen1 = arcLength(final_Points(1:halfIdx,1),final_Points(1:halfIdx,2));
arcLen2 = arcLength(final_Points(halfIdx:end,1),final_Points(halfIdx:end,2));

if arcLen1 < fullArcLen/2
   flg=1;
else
   flg=2;
end

k=0;
if flg==1
    while  arcLen1 < fullArcLen/2 && arcLen2 > fullArcLen/2
        arcLen1 = arcLength(final_Points(1:halfIdx+k,1),final_Points(1:halfIdx+k,2));
        arcLen2 = arcLength(final_Points(halfIdx+k:sz,1),final_Points(halfIdx+k:sz,2));
        k = k+1;
    end
    k = k-1;
    Index_crossing = halfIdx+k;
elseif  flg==2
      while  arcLen1 > fullArcLen/2 && arcLen2 < fullArcLen/2
        arcLen1 = arcLength(final_Points(1:halfIdx-k,1),final_Points(1:halfIdx-k,2));
        arcLen2 = arcLength(final_Points(halfIdx-k:sz,1),final_Points(halfIdx-k:sz,2));

        k = k+1;     
      end 
    k = k-1;
    Index_crossing = halfIdx-k;
end

quarterMark = round(Index_crossing/2); 
threequartersMark = round((Index_crossing + sz)/2);

middleGridPoints = zeros(3, 2, 'double');
middleGridPoints(1, 1:2) = final_Points(quarterMark, 1:2);
middleGridPoints(2, 1:2) = final_Points(Index_crossing, 1:2);
middleGridPoints(3, 1:2) = final_Points(threequartersMark, 1:2);
end