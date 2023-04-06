function themask = cleanMask(myimage, varargin)
% Input: a BW or logical mask that has too many small blobs, and so you
% would like to clean it up. This fucntion will return a new mask to you
% that has only blobs greater than a certain size. It also allows you to
% only include blobs of a certain minimum dsistance away from the largest
% blob. Preset values are minSize = 10,000 and minDist = 250. If the only
% input to the function is the mask, then these values of minSize and
% minDist will be used. 

% Michael Glendinning, 2023

minSize = 10000;
minDist = 250;
% if the function was called with additional values, they can overwrite the above
try minSize = varargin{1};
catch
end

try minDist = varargin{2};
catch    
end

regionData = regionprops('table', myimage,  {'Area', 'PixelList', 'PixelIdxList'});
regionData = sortrows(regionData,'Area','descend');

id = find(regionData.Area > minSize);
regionDataMat = regionData{id, 3}; %3 is the pixleList

idx = numel(id); 
for i = 1:idx
pix = regionDataMat{i}; % format is [y, x]
temp = polyshape(pix(:, 2), pix(:, 1), 'Simplify', true);
temp2 = temp.Vertices;
temp2(isnan(temp2(:, 1)), :) = []; %polyshape likes to add NaNs
myshapes{i} = temp2;
end

%make a list of every combination of pairings of id's to evaluate distances between

list1 = repmat([1:1:idx]', idx, 1); % (fyi transforms a matrix like [1,2,3,4] into : [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4])
list2 = repelem([1:1:idx]', idx); % (fyi transforms a matrix like: [1,2,3,4] into: [ 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4])
combos = horzcat(list1, list2); 

distances = zeros(size(combos, 1), 1);
for j = 1:size(combos, 1)
    if combos(j, 1) >= combos(j,2) 
        distances(j,1) = 0;
    else
        P = myshapes{combos(j, 1)};
        Q = myshapes{combos(j, 2)};
        D = pdist2(P, Q, 'euclidean', 'Smallest',1);
        distances(j,1) = min(D);
    end
end
distances = reshape(distances, [], idx);
distancesFromBiggest= distances(1, :);
id_keep = id(find(distancesFromBiggest < minDist));

themask = zeros(size(myimage, 1:2), 'logical');
for k = 1:numel(id_keep)
    idex = id_keep(k);
    list = regionData{idex, 2};
    themask(list{:}) = 1;
end

end

    
