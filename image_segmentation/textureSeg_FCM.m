function outOfBoundsMask = textureSeg_FCM(gray_img, colorImage, flag)
%outOfBoundsMask = textureSeg_FCM(gray_img, colorImage)
% Features:
% Peter Kovesi's log-gabor function for filtration
% LBP filter I implented.
% FCM algorithm fast and simple!
% Squared Hellinger Distance (works better than Euclidian)
filename = 'MSB_E7_1.mat';

if flag == 0
    rescaleFactor = 1200/max(size(colorImage, 1:2));
    colorImage = imresize(colorImage,rescaleFactor, {@oscResampling,4}); %gets resized after the switch-case block
    gray_img = imresize(gray_img,rescaleFactor, {@oscResampling,4}); %gets resized after the switch-case block
    imagesSave = struct('colorim', [], 'grayim', [], 'ID', []);
    
    imagesSave(1).colorim = colorImage;
    imagesSave(1).grayim = gray_img;
    imagesSave(1).ID = filename;
    
    save(filename, 'imagesSave', '-nocompression');
    
else
    imagesSave = load(filename);
    colorImage =   imagesSave.imagesSave(1).colorim;
    gray_img = imagesSave.imagesSave(1).grayim;
end
colorImage_lab = rgb2lab(colorImage);
gray2 = im2single(rescale(colorImage_lab(:,:,1), 0, 1));
gray2 = locallapfilt(gray2, 0.1, 0.5, 5);
gray2 = im2double(gray2);

% Apply local binary pattern to the image
lbp_im = local_binary_pattern(gray_img);
lbp_im2 = local_binary_pattern(gray2);

% Apply Gabor filter to the image
gabor_im = loggaborseg(gray_img);
gabor_im2 = loggaborseg(gray2);

% Ensure images are the same size before calling clustering FCN
sz = size(lbp_im, 1:2);
gabor_img = imresize(gabor_im, sz, 'bilinear');
gabor_img2 = imresize(gabor_im2, sz, 'bilinear');

% Combine the two filtered images using Fuzzy C-Means Clustering
[U2, ~] = FCM_clustering([lbp_im(:); gabor_img(:)], 2, 1.5);
[U3, ~] = FCM_clustering([lbp_im2(:); gabor_img2(:)], 2, 1.5);

[~,maxU2] = max(U2, [],1);
cluster1 = logical(reshape(maxU2==1, [], sz(2)));

[~,maxU3] = max(U3, [],1);
cluster2 = logical(reshape(maxU3==1, [], sz(2)));

montage({cluster1, cluster2});
outOfBoundsMask = logical(cluster1);

end

function lbp_im = local_binary_pattern(img)

% Set the LBP parameters
radius = 1;
num_neighbors = 8;
img = im2single(img);
% Preprocess the image
%    img = histeqfloat(img);
img = imresize(img, 0.8, 'bilinear');
img_padded = locallapfilt(padarray(img, [1, 1], 1, 'both'), 0.02, 1.75, 'NumIntensityLevels', 200);

% Define a function to calculate the LBP for a neighborhood
lbp_func = @(x) lbp(x, radius, num_neighbors);

% Calculate the LBP for each neighborhood
lbp_im = nlfilter(img_padded, [2*radius+1 2*radius+1], lbp_func);

% Remove the padding from the output image
lbp_im = lbp_im(radius+1:end-radius, radius+1:end-radius);

lbp_im = im2double(imresize(lbp_im, size(img, 1:2), 'bilinear'));

end

function gabor_im = loggaborseg(img_)

% Set the Log-Gabor filter parameters
theta = deg2rad(0:45:135);
scale = 3; %scale factor
orientation = pi/2; %orientation of the filter (in radians)
fo = 0.1; %centre frequency of the filter
BW = 0.3; %bandwidth of the filter

% Preprocess the image
img = img_; %imresize(img_, 1, 'bilinear');

% Create the Log-Gabor filter bank
gabor_bank = cell(length(theta), 1);
for p = 1:length(theta)
    gabor_bank{p} = loggaborfilter2(size(img), scale, orientation+theta(p), fo, BW);
    
end
gabor_bank = cellfun(@(x) imresize(x, size(img, 1:2) , 'bilinear'), gabor_bank, 'UniformOutput', false);

% Perform FFT on the image
F = fft2(img);

% Shift the zero-frequency component to the center of the array
Fsh = fftshift(F);

% Compute the filter response for each filter in the log-gabor filter bank
gabor_im = zeros([size(img), length(theta)]);
for p = 1:length(theta)
    % Compute the response in the frequency domain
    filter_fft = fft2(ifftshift(gabor_bank{p}));
    response_fft = Fsh .* filter_fft;
    response = ifft2(ifftshift(response_fft));
    % Compute the magnitude response
    gabor_im(:,:,p) = abs(response);
end

% Combine the magnitude responses of each filter to obtain the final filtered image
gabor_im = fftshift(gabor_im);
%gabor_im_sh = ifftshift(gabor_im);
% Convolve the filter bank with the image in the frequency domain
Fsh_filt = zeros([size(img), length(theta)]);
for p = 1:length(theta)
    % Compute the response in the frequency domain
    filter_fft = fft2(gabor_bank{p});
    response_fft = Fsh .* filter_fft;
    response = ifft2(ifftshift(response_fft));
    % Compute the magnitude response
    gabor_im(:,:,p) = abs(response);
    F_filt(:,:,p) = fftshift(Fsh_filt);

end
end

function lbp_code = lbp(neighborhood, radius, num_neighbors)
% Calculate the LBP code for a neighborhood

center_pixel = neighborhood(radius+1, radius+1);
lbp_code = 0;
for i = 1:num_neighbors
    % Calculate the coordinates of the neighbor pixel
    x = round(radius * cospi((i-1)*2/num_neighbors));
    y = round(radius * sinpi((i-1)*2/num_neighbors));
    neighbor_pixel = neighborhood(radius+1+y, radius+1+x);
    if neighbor_pixel >= center_pixel
        lbp_code = lbp_code + 2^(i-1);
    end
end

end

function [U_1, centroids_new, distMat] = FCM_clustering(intensityData, clusterN, fuzzifier_m)
% partition_matrix - updates the partition matrix and centroid values
% m is the fuzzifier parameter
% clusterN is the number of clusters
centroidP = regionprops('table', 'Centroid');
centroids_0 = centroidP.Centroid;

% Use k-means clustering to obtain initial centroids
[idx, centroids_0] = kmeans(intensityData, clusterN);

% Initialize U
Ucols = size(intensityData, 1);
U_0 = zeros(clusterN, Ucols);
for q = 1:Ucols
    [~, closestCentroid] = min(pdist2(centroids_0, intensityData(q,:), 'euclidean'));
    U_0(closestCentroid, q) = 1;
end


n=1;
while n<16
    % Calculate new centroids
    u_m = U_0 .^ fuzzifier_m;
    centroids_new = (u_m * intensityData) ./ sum(u_m, 2);
    
    % Calculate distance matrix using the Squared Hellinger Distance
    distMat = pdist2(centroids_new, intensityData, @(x, y) 0.5*sum((sqrt(x) + sqrt(y)).^2));
    
    % Remove zero distances from same_centroid
    same_centroid = find(distMat == 0);
    [same_centroid_idx, same_centroid_cluster] = ind2sub(size(distMat), same_centroid);
    distMat(same_centroid) = NaN;
    
    % Calculate reciprocal of distance matrix
    recip_distMat = 1 ./ distMat;
    recip_distMat(same_centroid) = 0;
    
    % Calculate partition matrix
    rD = recip_distMat.^(fuzzifier_m - 1);
    U_1 = rD ./ sum(rD, 1);
    U_1(:, same_centroid_cluster) = 0;
    for i = 1:length(same_centroid_idx)
        U_1(same_centroid_idx(i), same_centroid_cluster(i)) = 1;
    end
    
    % Check for convergence
    if any(abs(U_0 - U_1) > 0.005)
        n=n+1;
        U_0 = U_1;
    else
        n=16;
    end
end


end


function f = loggaborfilter2(N,scale,orientation,fo,phaseoffset)
%Peter Kovesi

[X, Y] = meshgrid(-fix(N/2):fix(N/2), fix(N/2):-1:fix(-N/2));

% Radial filter component
Xp = X * cos(orientation) + Y * sin(orientation);
Yp = -X * sin(orientation) + Y * cos(orientation);
r = hypot(Xp, Yp);

% Log-Gabor filter
f = exp( (-(log(r/fo)).^2) / (2 * log(scale)^2) );
f = f .* cospi(2 * (Xp/fo) + phaseoffset);

% Set the DC component of the filter to zero (Gaussian high-pass)
f(fix(N/2)+1,fix(N/2)+1) = 0;

end

