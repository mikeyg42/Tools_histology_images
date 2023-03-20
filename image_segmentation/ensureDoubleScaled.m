function output = ensureDoubleScaled(varargin)
% usage: [output] = ensureDoubleScaled(im_1, (TF))
% OR: {output} = ensureDoubleScaled(im_1, im_2,..,im_N, (TF))

% you can list one or multiple images... listing multiple images is a bit 
% faster than repeatedly calling fcn. If you do input multiple images,
% output will be a cell array of images, such that im_k--> output{k}. but, if
% only one image is in input, output will not be a cell array, instead it will
% simply be an image.

% OPTIONAL FLAG: there is option to inlude as an input a logical value, TF.
% -- TF can be omitted, TRUE, or FALSE.
% -- it MUST go after all the images.
% -- when =TRUE, the function will use RESCALE to ensure that the range of 
% intensity values is at least 0.6 (i.e. so that maxValue-minValue>=0.6). 
% -- You can omit TF, it will default to FALSE.

% Procedure:
% 1. converts to floating with im2double. 
% 2. Then it rescales values to be within range of [0,1].
% 3. it checks for any complex values. If we one is found, 
% first we attempt to linearly interpolate with a neighboring pixels. 
% In absense of a real valued neighboring pixel, we can relace with the 
% median of the absolute values of its 7x7 local neighborhood. 

% NOTE: NaNs and Inf's are left AS IS. Address seperately if necessary.

% Michael Glendinning, 2022

%% step 1: parse input
if islogical(varargin{nargin})
    if varargin{nargin}   %ie if true
        flag = 1;
    else
        flag = 0;
    end
    numIms=nargin-1;
    myIms = varargin(1:numIms);
else 
    flag = 0;
    numIms = nargin;
    myIms = varargin;
end

if nargout ~=1 
    error('ensureDoubleScaled: change output of ensureDoubleScaled to be only 1 variable');
end

%% step 2: begin looping through input images
if numIms >1
    output = cell(1,numIms);
end

for qq = 1:numIms
    
    img = myIms{qq};
    if ~isnumeric(img)
        error('ensureDoubleScaled: its required that every input be an image, ie a numerical array, (except for the final optional logical input)')
    end
    
%% step 3: convert to double using IM2DOUBLE    
    if ~isa(img, 'float')
        img = im2double(img);
    end

%% step 4: ensure that range of values is between 0 and 1    
% you need these steps because if your input is a float, im2double
 % won't rescale, but if they are out of range you need to fix!
 
    minVal = min(img(:));
    maxVal = max(img(:));
    
   if maxVal>1.1
       rescale(img, minVal, 1);
       maxVal = 1;
   end
   
   if minVal<0
       if (1-maxVal) > abs(minVal)
           img = imadd(img,double(abs(minVal)));
       else
           img = rescale(img,0, maxVal);
       end
       minVal = 0;
   end
   
   % WHEN INDICATED BY TF: Ensure that range(ie max-min) of image is >/= 0.6 
   if maxVal-minVal<0.6 && flag == 1
       disp('maxVal-minVal<0.6 and so I increased contrast slightly');
       avgVal = (maxVal+minVal)/2;
       %here I ensure that the intensty range is always at least 0.6 
       if avgVal<0.7 && avgVal>0.3
           img = rescale(img, avgVal-0.3,avgVal+0.3);
       elseif avgVal>0.7
           img = rescale(img, 0.4,1);
       elseif avg<0.3
           img = rescale(img, 0,0.6);
       end
   end
    	
% very infrequently, I will suffer from a few pixels randomly being complex
% values? this is a rather inelegant way to fix that problem

    if any(imag(sum(img))~= 0) %not entirely sure why, but its 25x as fast using sum... (imag(img(:)) = 0.49sec vs imag(sum(img)) = 0.018sec)
        img = convertImaginaryPixels2Real(img);
    end

    % populate output variable called OUTPUT 
    if numIms > 1
    output{qq} = img;
    elseif numIms == 1 
    output = img;
    end
end


end

function img = convertImaginaryPixels2Real(img)
   warning('ImageError:ComplexNumberIssue', ...
'Array contains complex numbers... Converting pixels using nearby values')

    if ndims(img)==3
        sz = size(img);
        img = reshape(img, sz(1), [], 1);
        threeD_flag = 1;
    else
        threeD_flag = 0;
    end

    [rows, cols, vals] = find(imag(img));
    imgi = padarray(img, 1, 1+1i, 'both'); %pad with 1+i so that you can easily find i 4/8 connective neighborhood
    rows = rows+1; cols = cols+1; %adjust for the padding

    for pp = 1:length(vals)
        r = rows(pp); c= cols(pp);
        nhood = [imgi(r-1,c); imgi(r+1, c); imgi(r,c-1); imgi(r, c+1)];
        replacementVals = find(nhood==real(nhood)); %get the indexes for the nonnegative values in 4 point conn nhood
        if ~isempty(replacementVals)
            imgi(r,c) = sum(nhood(replacementVals))/numel(replacementVals);
        else
            nhood2 = [imgi(r-1,c-1); imgi(r+1, c+1); imgi(r+1,c-1); imgi(r-1, c+1)]; %extend to 8point connnectivity
            replaceVal2= find(nhood2==real(nhood2));
            if ~isempty(replaceVal2)
            imgi(r,c) = sum(nhood(replaceVal2))/numel(replaceVal2);
            else
                %if you get to this point, there are so many complex values that it seems essential to consider the 
                medval = median(abs(imgi(r-3:r+3, c-3:c+3)), 'all'); %evaluate a 6x6 neighborhood median(abs(nhood)))
                imgi(r,c) = min(max(0, medval), 1); %ensure its less than 1 and greater than 0
            end
        end
    end
    %remove padding
    img = imgi(2:end-1, 2:end-1);

    if threeD_flag == 1
        img = reshape(img, sz);
    end
end
