function imAdjRGB = preprocessRawRGBims(myImage, scaleFactor)

%sometimes a 4th empty channel is added by Zen Blue (an alpha channel I think)... a weird Zeiss software glitch
sz = size(myImage, 1:3);
if sz(3)>3 
    myImage(:,:,4:end)=[];
end

% this converts to double and rescaled everything do fall in [0,1] range.
myImage = ensureDoubleScaled(myImage, true);

%downsampling image by scale factor. don't scale down too much!! just
%enough that code will run....
sz = size(myImage, 1:2);
newSize = ceil([sz(1), sz(2)].*scaleFactor);
myImage = imresize(myImage, newSize,{@oscResampling, 4}); %downsampling

imAdjRGB = backgroundCorrectRGB(myImage, 'RGB_stacked'); %removing background and mosaic artifacts
end