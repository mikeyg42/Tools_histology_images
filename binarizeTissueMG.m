function BinMask = binarizeTissueMG(RGBim, varargin)

%% parsing optional input of prespecified polybuffer amount (
flag = 1;
try flag = isnumeric(varargin{1});
catch
    flag = 0;
end
if flag == 1
    bufferAmount = varargin{1};
else
    bufferAmount = 9;
end

%% relies on imfill, imdiffusefilt, and the roicolor on the luminosity chl (L*A*B image)

lab = rgb2lab(RGBim);
lscaled = rescale(lab(:,:,1), 0, 1);
lCh = imcomplement(lscaled);

se = strel('disk', 21, 8);
imGray = imfill(imclose(lCh, se));
% this imclose and fill helps bring the lesion (much less luminance) up to the threshold for binarization, at the cost of some precision

%basic binarization of gray image
[Tlevel, ~] = graythresh(imGray);
Idiffusion = imdiffusefilt(imGray);
newROImask = roicolor(Idiffusion, (Tlevel/1.5), 1);

% select the largest blob
BW = bwareafilt(logical(newROImask), 1);

%define the boundaries using contour
ct = contourc(BW, [0.5, 0.5]);
conX = ct(1,2:end);
conY = ct(2,2:end);

%convert contour into a polyshape and simplify the edges
ps = polyshape(conX, conY, 'Simplify', true);

%buffer the edges a bit to smoothen some jagged edges caused by the fat SE.
%remake your mask
ps2 = polybuffer(ps.Vertices, 'points', bufferAmount);
bufferedMask = poly2mask(ps2.Vertices(:, 1), ps2.Vertices(:, 2), size(RGBim, 1), size(RGBim, 2));

target = imcomplement(imgradient(lscaled));
%remove your padding of edges with 30 iterations of active contour
BinMask = activecontour(target, bufferedMask, 30);
end

