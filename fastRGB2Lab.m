function [L,a,b] = fastRGB2Lab(rgbIm, disp_y_n)
% RGB2LAB Convert an image from RGB to CIELAB

% disp_y_n = 1 --> figure will appear showing each channel
% disp_y_n = 0 --> skip the visualization

% RGB2Lab takes a single M x N x 3 image, 
% and returns an image in the CIELAB color space individually.  
% RGB values can be double or uint8.
% Values for L are in the range [0,100] (double)
% Values for a and b are roughly in the range [-110,110] (double)
%
% This transform is based on ITU-R Recommendation BT.709 using the D65
% white point reference. (MATLAB native imp uses D50??) 
% The error in transforming RGB -> Lab -> RGB is approximately 10^-5.  
%
% By Mark Ruzon from C code by Yossi Rubner, 23 September 1997.
% (their last update -- 30 March 2009. forked from file exchange)
% I've tweaked it slightly here for my needs (MG - Dec 2022)

% I've found that particularly close to 0 the L channel is quite off 
% at least compared to MATLABs function. a rescale(L, 0,100) fixes, at cost
% of some speed of course.

[R, G, B] = imsplit(rgbIm);
RGB = [R(:)'; G(:)'; B(:)'];

RGB = double(RGB);
if any(RGB(:) > 1.0)
  RGB = RGB./255;
end

% RGB to XYZ
MAT = [0.412453 0.357580 0.180423;
       0.21267 0.715160 0.072169;
       0.019334 0.119193 0.950227];
XYZ = MAT*RGB;

% Normalize for D65 white point
X = XYZ(1,:)./0.950456;
Y = XYZ(2,:);
Z = XYZ(3,:)./1.088754;

% Set a threshold
T = 0.008856;

XT = X > T;
YT = Y > T;
ZT = Z > T;

Y3 = Y.^(1/3); 

fX = XT .* X.^(1/3) + (~XT) .* (7.787 .* X + 16/116);
fY = YT .* Y3 + (~YT) .* (7.787 .* Y + 16/116);
fZ = ZT .* Z.^(1/3) + (~ZT) .* (7.787 .* Z + 16/116);

[m, n] = size(rgbIm, 1:2);
L = reshape(YT .* (116 * Y3 - 16.0) + (~YT) .* (903.3 * Y), m, n);
a = reshape(500 * (fX - fY), m, n);
b = reshape(200 * (fY - fZ), m, n);

if disp_y_n == 1
    figure;
    subplot(1, 3, 1);
    imshow(L, [], 'Border', 'tight');
    title('L Image');
    subplot(1, 3, 2);
    imshow(a, [], 'Border', 'tight');
    title('A Image');
    subplot(1, 3, 3);
    imshow(b, [], 'Border', 'tight');
    title('B Image');
end
