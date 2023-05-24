function symD = imFractionDeriv(img,v)
% symD = imFractionDeriv(img, v)
% Computes the left and right Riemann-Liouville fractional derivatives of an image
% using the Caputo-Grümwald method with the order of differentiation v. Note that this is  
% neither the Grumwald-Letnikov derivative nor the Caputo derivative. Rather, it is the
% numerical method of calculating the Riemann-Liouville derivative described by Caputo and 
% Grumwald which increases the accuracy/efficiency of the numerical calc as initially described.
% The methods produces a LEFT and a RIGHT derivative, very similar to the the forward vs.
% the backwards difference re: integer derivatives, i.e. each describes the behavior of curve 
% on one half or the other of a point. If one has a application that is concerned much more 
% with a particular side of each point, then perhaps theyd be well-served to not calculate the 
% arithmetic mean of the L and R derivatives (which I've done to find a balanced emphasis on 
% either sides of a point) as I do in the last line of this code. 
% 
% Note that v must be positive. It will give valid result for v = 0, which is a special a case 
% where function is reduced to the dirac delta function. v is recommended not to exceed 2 or 3,
% which would be considered the general ceiling for approximations of higher-order
% fractional derivatives... this method is woefully inaccurate with anything much above 2.
%
% Note also that img MUST be floating point precision or an error will be thrown. 

% Michael Glendinning, 2023
% ======================================================================== %

assert(isfloat(img)&& ismatrix(img));
assert(v>=0)

[rws,cols] = size(img);

% Create grid of x and y coordinates
[xgrid,ygrid] = meshgrid(1:cols,1:rws);

% Compute L1 distance between neighboring pixel values 
% (for images with a pixel aspect ratio = 1, this will always be
% an array full of 2's, with 0's on the last row/col. 
% Thus, always it'll be that: N = 2 and h =1.99ish
distanceMap = abs(xgrid - circshift(xgrid,[0,-1])) + abs(ygrid - circshift(ygrid,[-1,0])); % L1 distance map

% last row+col inherently must always be zero...
distanceMap(rws,:) = 0; 
distanceMap(:,cols) = 0; 

% Compute weights for Caputo-Grümwald method (using a power-law kernel)
h = mean(distanceMap(:));
N = ceil(2/h); % for most pictures (pix aspect ratio = 1), N will always be 2
weights = zeros(1,N+1);

% "gamma" is the special function described by Euler, namely an Euler integral of the 2nd kind) 
% when x>0 its given by gamma(n+1) = n! , and notable gamma(0.5) = pi^0.5
for k = 0:N
    weights(k+1) = (-1)^k * gamma(v+1) / (gamma(k+1) * gamma(v-k+1)); 
end

% Compute L and R Riemann-Liouville fractional derivatives by looping column-wise through img
D_left = zeros(rws,cols); %preallocate
D_right = zeros(rws,cols);

repWeights = repmat(weights, rws, 1); %pre-calculate this!
for g = 1:cols
    if g <= N
        D_left(:,g) = sum(repmat(weights(1:g+1), rws, 1) .* img(:,1:g+1) .* distanceMap(:,1:g+1).^(v+1),2);
        D_right(:,g) = sum(repWeights .* img(:,g:g+N) .* distanceMap(:,g:g+N).^(v+1),2);
        
    elseif g > cols-N
        D_left(:,g) = sum(repWeights .* img(:,g-N:g) .* distanceMap(:,g-N:g).^(v+1),2); 
        D_right(:,g) = sum(repmat(weights(1:cols-g+1), rws, 1) .* img(:,g:cols) .* distanceMap(:,g:cols).^(v+1),2);
    
    else 
         D_left(:,g) = sum(repWeights.* img(:,g-N:g) .* distanceMap(:,g-N:g).^(v+1), 2);
        D_right(:,g) = sum(repWeights .* img(:,g:g+N) .* distanceMap(:,g:g+N).^(v+1),2);
    end
end
D_left = D_left / h^(v+1);
D_right = D_right / h^(v+1);

% Take the arithmetic mean of the left and right Riemann-Liouville fractional derivatives
symD = (D_left+D_right)./2;

end