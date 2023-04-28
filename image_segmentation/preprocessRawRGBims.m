function imAdjRGB = preprocessRawRGBims(myImage, scaleFactor)
% imAdjRGB = preprocessRawRGBims(myImage, scaleFactorToResize)

%sometimes a 4th empty channel is added by Zen Blue (an alpha channel I think)... a weird Zeiss software glitch
sz = size(myImage, 1:3);
if sz(3)>3 
    myImage(:,:,4:end)=[];
end

% this converts to double and rescaled everything do fall in [0,1] range.
myImage = ensureDoubleScaled(myImage, true);

%downsampling image by scale factor. shink too much!! just enough that code will run....
if scaleFactor ~=0
newSize = ceil([sz(1), sz(2)].*scaleFactor);
myImage = imresize(myImage, newSize,{@oscResampling, 4}); %downsampling
end

imAdjRGB = backgroundCorrectRGB(myImage, 'RGB_stacked'); %removing background and mosaic artifacts
end