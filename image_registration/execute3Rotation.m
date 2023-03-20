function [rotatedImage, rotatedMask, rotCorner] = execute3Rotation(myImage, myImageMask, myCorners, theta)
%[rotatedImage, rotatedMask, rotCorner] = execute3Rotation(myImage, myImageMask, myCorners, theta)
%rotate all at once theta degrees (not radians) an image, its mask, and the
% 4x2 matrix of corner points. 

% Images are translated so that their centers are on the origin. Then they
% are rotated, using spline interpolation. Finally, they are rotated back
% to their original location. 

% Michael Glendinning, 2023

if ~isfloat(myImage)
myImage = im2double(myImage); end

if ~isfloat(myImageMask)
myMask = im2double(myImageMask); end

rws = size(myImage, 1);
cols = size(myImage, 2);
[xi,yi] = meshgrid(1:cols,1:rws);  

xOriginTranslate = cols/2;
yOriginTranslate = rws/2;

deltaX = xi - xOriginTranslate;
deltaY = yi - yOriginTranslate;

x_interpolated =deltaX.*cosd(theta) + deltaY.*sind(theta) + xOriginTranslate;   % transformed coordinates (new pixel locations)
y_interpolated =deltaY.*cosd(theta) - deltaX.*sind(theta) + yOriginTranslate;

rotatedImage = interp2(xi, yi, myImage, x_interpolated, y_interpolated, 'spline');  
rotatedMask = interp2(xi, yi, myMask, x_interpolated, y_interpolated, 'spline');

% I've opted below for pre-mulitply convention for the points
% instead of post-multiply as I did above, please don't be confused.

%translate to the origin
translateInv = eye(3);
translateInv(1:2, 3) = [-xOriginTranslate; -yOriginTranslate];
dataMatrix = [myCorners, ones(length(myCorners), 1)];
myCorn = translateInv*dataMatrix';
myCorn(3, :)= [];

%rotate about the origin
newCorn = zeros(size(myCorn));
for j = 1:length(myCorn)
newCorn(1:2, j) = [cosd(theta), -sind(theta); sind(theta), cosd(theta)]*myCorn(1:2, j);
end
rotCorn = newCorn';

%translate back to the original location
dataMatrix2 = [rotCorn, ones(length(rotCorn), 1)];
translateFwd = eye(3);
translateFwd(1:2, 3) = [xOriginTranslate; yOriginTranslate];
rotCo = translateFwd*dataMatrix2';
rotCo(3,:)= [];
rotCorner = rotCo';

% C1 = imfuse(rotatedImage,myImage,'falsecolor','Scaling','none','ColorChannels',[1 2 0]);
% C2 = imfuse(rotatedMask,myImageMask,'falsecolor','Scaling','none','ColorChannels',[1 2 2]);
% C3 = imfuse(rotCorner,myCornerPointVisual,'falsecolor','Scaling','none','ColorChannels',[2 1 2]);
% figure; montage({C1,C2,C3})
end




