function outputRGBim = evenFasterCorrection(RGBfloat)

se = strel('disk', 12, 8);

assert(isfloat(RGBfloat));
RGBpad = padarray(RGBfloat, [3, 3], 'replicate', 'both');
    
smallPADsz = size(RGBpad, 1:3);
onesMat = repmat(1, smallPADsz(1), smallPADsz(2)*3); % repmat is twice as fast as ones([m,n])
    
inverted2Dim_unfilt = onesMat - cat(2, RGBpad(:,:,1) , RGBpad(:,:,2) , RGBpad(:,:,3)); %dim=2, like a hot dog, not a burger
    
invFiltered2D = inverted2Dim_unfilt -  imopen(inverted2Dim_unfilt, se);
filtered2Dim =  onesMat - invFiltered2D; 

ncol = smallPADsz(2);
RGBadj_pad = zeros(smallPADsz, 'like',RGBfloat);
RGBadj_pad(:,:,1) = filtered2Dim(:, 1:ncol);
RGBadj_pad(:,:,2) = filtered2Dim(:, 1+ncol:2*ncol);
RGBadj_pad(:,:,3) = filtered2Dim(:, 1+2*ncol:end);

outputRGBim = RGBadj_pad(4:end-3, 4:end-3, :);

end

