function imgFin = mikeMedianFilter(img_, numIter, maxSizeDimension, colorspace)
%syntax: imgFin = mikeMedianFilter(img_, numIter, maxSizeDimension, colorspace)
% numIter = 1 usually! in order to make multiple iterations more effective,
% I've included 3 slightly different 5x5 kernels that are cycled between
% each iteration.

% maxSizeDimension is an integer < 2500 (otherwise SLOW! 1800) works very
% fast. Image will be downsampled to match this value, then upsampled back
% to original after. To avoid any resizing, input 0 for maxSizeDimension.

% colorspace == {'RGB', or 'Luminance' or 'grayscale'}... for Luminance,
% the input should stilll be RGB. Will convert to RGB, median filter the first
% channel, then concatenate with the unchanged channels 2 and 3! if your image 
% is already converted luminance, indicate grayscale! 

% IMG_ should be a DOUBLE always. if not single or double float, 
% it will be converted to double before doing median filter. it resizes back 
% up to OG size.  
% weighted averages is taken with the average image in order to reduce any
% artifacts from resampling or filtering

if ~isfloat(img_)
    img_ = ensureDoubleScaled(img_);
end
img_ = squeeze(img_);
colorspace = upper(colorspace(1)); %either R or L or G

if ndims(img_) < 3
    if ~strcmp(colorspace,'G')
    warning('Colorspace is incorrectly defined.  Assuming input is a grayscale image because it has only 2 dimensions');
    end
    img = cat(3, img_, img_, img_);
elseif ndims(squeeze(img_)) == 3
    img = img_;
else
    error('must be RGB image! too many dimensions in image for median filter');
end

[nrow, ncol] = size(img,1:2);
if nrow>ncol %we want there to be more columns than rows! 
    img = rot90(img);
    [nrow, ncol] = size(img,1:2);
    flag = 1;
else
    flag = 0; 
end

if maxSizeDimension == 0
   outIm_row = nrow;
   outIm_col = ncol;
else
    outIm_col = maxSizeDimension; % value from input of function
    outIm_row = ceil(nrow*outIm_col/ncol);
end

% pre-allocate 
imgFin = zeros(outIm_row,outIm_col,3, class(img));

% shrink image to hasten the very slow median filter
imgResized = imresize(img, [outIm_row,outIm_col] ,{@oscResampling,4}); 

% median filter to blur
filt1 =  [0 1 1 1 0; 
         1 1 1 1 1; 
         1 1 1 1 1; 
         1 1 1 1 1; 
         0 1 1 1 0];
filt2 = ones(5); %after the first iteration it iwll just continue 
filt3 = [0 0 1 0 0; 
         0 1 1 1 0; 
         1 1 1 1 1; 
         0 1 1 1 0; 
         0 0 1 0 0];
     
     switch colorspace
         
         case 'R' %RGB space
             
             [a1 , a2 , a3] = imsplit(imgResized);
             
             jj = 0;
             while jj ~= numIter
                 if mod(jj, 3) ==0, filt = filt1; num=11; 
                 elseif mod(jj,3)==1, filt= filt2; num = 13; 
                 else, filt = filt3; num = 7; 
                 end
                 
                 a1 = ordfilt2(a1, num, filt, 'symmetric'); % can brighten image by choosing 12 or higher
                 a2 = ordfilt2(a2, num, filt, 'symmetric');
                 a3 = ordfilt2(a3, num, filt, 'symmetric');
                 jj = jj+1;
             end
             
             imgFin(:,:,1) = a1;
             imgFin(:,:,2) = a2;
             imgFin(:,:,3) = a3;
             
         case 'L' %Luminance
             
             [a1, a2, a3] = imsplit(rgb2lab(imgResized));
             
             a1 = a1./100;
             jj = 0;
             while jj ~= numIter
                 if mod(jj, 3) ==0, filt = filt1; num=11; 
                 elseif mod(jj,3)==1, filt= filt2; num = 13; 
                 else, filt = filt3; num = 7; 
                 end
                 jj = jj+1;
                 a1 = ordfilt2(a1, num, filt, 'symmetric');
             end
             L2 = a1.*100;
             L3 = cat(3, L2, a2, a3);
             
             imgFin = lab2rgb(L3);
             
         case 'G' %grayscale
             try a1 = rgb2gray(imgResized);
             catch
                 if 2 == ndims(imgResized)
                     a1 = imgResized;
                 else 
                     disp(strcat(class(imgResized), ' Is unable to be converted grayscale'));
                 end
             end
                
             jj = 0;
             while jj ~= numIter
                 if mod(jj, 3) ==0
                     filt = filt1; num=11; 
                 elseif mod(jj,3)==1
                     filt= filt2; num = 13 ;
                 else
                     filt = filt3; num = 7;
                 end
                 jj = jj+1;
                 a1 = ordfilt2(a1, num, filt, 'symmetric'); % can simultaneously brighten image by choosing num values that are slightly higher...
             end
             imgFin(:,:,1) = a1;
             imgFin(:,:,2:3) = []; 
     end
     
     if flag == 1
         imgFin = rot90(imgFin, -1);
         img = rot90(img, -1);
     end
     
     imgFin = imresize(imgFin, size(img_,1:2), {@oscResampling,4});
     
     imgFin = (img+imgFin.*2)./3; %if this doesn't work when tried, its because of the dimension mismatch from rot90
     
end

