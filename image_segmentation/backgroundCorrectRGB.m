function RGBadj = backgroundCorrectRGB(inputImage, descriptor)
% adjustedRGB = backgroundCorrectRGB(inputImage, descriptor)
% 
% "descriptor" must be a string saying: 
% 1. 'Gray',        2. 'RGB_stacked',      3. 'RGB_tiled'. 

% - RGB_Stacked refers to normal RGB image with 3 pages.  
% - RGB_Tiled indicates that the input image was made into a montage
%    of the 3 channels adjacent eachother in 2D space. 
% - Stacked images will be converted to tiled images as the first step 

% If nothing is included in "descriptor", it will default to RGB_stacked
% for 3D arrays, and "Gray" for 2D arrays

% Nonuniform background is corrected by the classic morphology manipulation
% of subtracting an "opened image" from the original. This manipulation 
% is also called the tophat filter.

% Because of the morphology functionality requiring a black background with 
% white foreground, we take the complement before and after invoking the tophat function.
%
%           As such, the general formula is:
% bckgndGOOD = imcomplement(IMTOPHAT(imcomplement(bckgndBAD),se))
% 
% I've also included a padding step to the RGB images. I'm not sure its
% necessary though, but it cannot hurt! 
%
%
% Michael Glendinning, 2022.
%======================================================================

% regardless of the input image type, we use this structuring element
% for our tophat filtering 
se = strel('disk', 12, 8);

if nargin == 1
    if ndims(inputImage) == 3
    descriptor = 'RGB_stacked';
    else
    descriptor = 'Gray';
    end
elseif ~isstring(descriptor) && ~ischar(descriptor)
    disp('format of input named DESCRIPTOR is wrong. must be char or string');
    return
end
      
if ~isa(inputImage,'float')
    RGB = ensureDoubleScaled(inputImage);
else 
    RGB = inputImage;
end

if strcmpi(descriptor, 'Gray')
% it facilitates the tophatfilter to have more contrast (thus imadjust).
    inverted2Dim_unfilt = 1 - imadjust(RGB); %... however we won't keep it that way...
    filtered2Dgray = 1 - fastTophat(inverted2Dim_unfilt,se);
    RGBadj = rescale(filtered2Dgray, min(RGB(:)), max(RGB(:)));
     %... we rescale back to original contrast here ^^.
 
else


if strcmpi(descriptor, 'RGB_stacked')
    flag = 0;%flag reminds us not to forget to remove PAD 
    
    RGBpad = padarray(RGB, [3, 3], 'replicate', 'both');
    
    small_sz = size(RGBpad, 1:3);
    onesMat = repmat(1, small_sz(1), small_sz(2)*3);% twice as fast as ones([m,n])
    
    inverted2Dim_unfilt = onesMat - cat(2, RGBpad(:,:,1) , RGBpad(:,:,2) , RGBpad(:,:,3)); %dim=2, like a hot dog, not a burger
    
elseif strcmpi(descriptor, 'RGB_tiled')
    flag = 1; %flag informs us whether or not there is a padding to remove. 1 = NO pad
    
    big_sz = size(RGB, 1:2);
    small_sz = [big_sz(1), big_sz(2)/3, 3];
    
    onesMat = repmat(1, small_sz(1), small_sz(2)*3);
    inverted2Dim_unfilt = onesMat-RGB; %inverse

else
    disp('incorrect 2nd input. 2nd input should be a string describing input: Gray, RGB_stack, or RGB_tiled');
    return
end
fastTH = fastTophat(inverted2Dim_unfilt, se);
filtered2Dim =  onesMat - fastTH; 

ncol = small_sz(2);
RGBadj_pad = zeros(small_sz, 'double');
RGBadj_pad(:,:,1) = filtered2Dim(:, 1:ncol);
RGBadj_pad(:,:,2) = filtered2Dim(:, 1+ncol:2*ncol);
RGBadj_pad(:,:,3) = filtered2Dim(:, 1+2*ncol:end);

if flag == 0
RGBadj = RGBadj_pad(4:small_sz(1)-3, 4:small_sz(2)-3, 1:3);
else
    RGBadj = RGBadj_pad;
end
end
end

function img2 = fastTophat(img, SE)
% this is the "white top-hat filter" 
img2 = img -  imopen(img,SE);
end
