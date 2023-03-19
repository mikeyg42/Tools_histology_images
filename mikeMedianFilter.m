function imgFin = mikeMedianFilter(img_, numIter, maxSizeDimension, colorspace)
% numIter= 1 usually!
% maxSizeDimension < 2500 otherwise SLOW! 1800 works very fast. 
% colorspace == {'RGB', or 'Luminance' or 'grayscale'};
% IMAGE_ should be a DOUBLE always
% after doing median filter. it resizes back up to OG size and then
% weighted averages back in the original image. 

img = ensureDoubleScaled(img_);
colorspace = upper(colorspace(1)); %either R or L or G
if ndims(img) < 3
    if ~strcmp(colorspace,'G')
    warning('Colorspace is incorrectly defined.  Assuming input 2-dimensional array is a grayscale image');
    end
    img = cat(3, img, img, img);
end

flag = 0;
[nrow, ncol] = size(img,1:2);
if nrow>ncol %we want there to be more columns than rows! 
    img = rot90(img);
    [nrow, ncol] = size(img,1:2);
    flag = 1;
end

outIm_col = maxSizeDimension;
shrinkFactor = outIm_col/ncol;
outIm_row = ceil(nrow*shrinkFactor);

% pre-allocate 
imgResized = zeros(outIm_row,outIm_col,3, 'double');
imgFin = zeros(outIm_row,outIm_col,3, 'double');

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
                 if mod(jj, 3) ==0, filt = filt1; num=11; elseif mod(jj,3)==1, filt= filt2; num = 13; else, filt = filt3; num = 7; 
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
             
             [a1, a2, a3] = fastRGB2Lab(imgResized, 0);
             
             a1 = a1./100;
             jj = 0;
             while jj ~= numIter
                 if mod(jj, 3) ==0, filt = filt1; num=11; elseif mod(jj,3)==1, filt= filt2; num = 13; else, filt = filt3; num = 7; 
                 end
                 jj = jj+1;
                 a1 = ordfilt2(a1, num, filt, 'symmetric');
             end
             L2 = a1.*100;
             L3 = cat(3, L2, a2, a3);
             
             imgFin = lab2rgb(L3);
             
         case 'G' %grayscale
             
             [a1 , ~ , ~] = imsplit(imgResized);
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
                 a1 = ordfilt2(a1, num, filt, 'symmetric'); % can brighten image by choosing 12 or higher  
             end
             imgFin(:,:,1) = a1;
             imgFin(:,:,2:3) = [];
             img(:,:,2:3) = [];
     end
     
     imgFin = imresize(imgFin, size(img,1:2), {@oscResampling,4});
     
     FlagPeep = 0;
     try imgFin = (img+imgFin.*2)./3; %if this doesn't work when tried, its definitely because of the dimension mismatch from rot90
     catch
         FlagPeep = 1;
         imgFin = rot90(imgFin, -1);
         imgFin = (img+imgFin.*2)./3;
     end
     
     if FlagPeep == 0 && flag == 1
         imgFin = rot90(imgFin, -1);
     end

     
end

