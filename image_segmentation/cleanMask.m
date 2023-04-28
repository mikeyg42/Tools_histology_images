function myCleanMask = cleanMask(mymask, varargin)
% syntax1: myCleanMask = cleanMask(myimage, minSize, minDist);
% syntax2: myCleanMask = cleanMask(myimage); 
% syntax3: myCleanMask = cleanMask(myimage, NaN, minDist);

% Goal: simple method of cleaning-up the output of any "crude" or imperfect 
% segmentation algorithm, treating each connected component of the binary input image
% as a blob; blobs are evaluated based on their area and relative distance to the largest
% blob. 
% Inputs: a BW or logical mask that potentially has too many smaller blobs, and you
% want to potentially clean it up. 
% Output: new mask,i.e. 'myCleanMask', that still contains all the blobs w/ 
% areas > minSize, but those blobs smaller than minSize are gone. Also, the blobs 
% that were not within a certain distance (the minDist) relative to the largest blob,
% are removed. 
% ====================================
% Pre-set values are minSize = 10,000 and minDist = 250. If the only
% input to the function is the mask, then these defaults for minSize and
% minDist will be used. To specify one and not the other, use NaN; any NaNs in varargin
% will be replaced with default values.

% Michael Glendinning, 2023

% Parse inputs:
if nargin == 1
    minDist = 250;
    minSize = 10000;
elseif nargin > 1
    if ismissing(varargin{1}) || isempty(varargin{1})
        minDist = 250;
    elseif isnumeric(varargin{1})
        minDist = varargin{1};
    else
        error('incorect inputs to clearMask function varargin{1}');
    end
    
    if nargin > 2    
        if ismissing(varargin{2}) % this would be the weird instance of writing clearMask(mask, 9, NaN);
            minSize = 10000;
        elseif isnumeric(varargin{2}) 
            minSize = varargin{2};
        else
            error('incorect inputs to clearMask function varargin{2}');
        end
    elseif nargin == 2
         minSize = 10000;
    end
end

% evaluate the largest blob
largest1blobIm = bwareafilt(mymask, 1);
area = sum(largest1blobIm, 'all');

% If the largest blob is smaller than minSize value, final result will be "overcleaned"!!
if minSize > area
    myCleanMask = mymask;
    disp('could not clean mask because no blob was larger than minSize... mask returned as is');
    return;
end
% Filter all blobs except those with areas between minSize and +infinity 
mymask = bwareafilt(mymask, [minSize, Inf]);

% dilate the largest blob with a structural element the size of minDist
SE_radius = strel('disk', minDist,8) ;
largestblob_buffered = imdilate(largest1blobIm, SE_radius);

% Any blob that intersects with the dilated largest blob we want to keep. We can recover
% each of those blobs with reconstruction, using the intersection as our seed. 
inRangeBlob_seed = largestblob_buffered & mymask;
myCleanMask = imreconstruct(inRangeBlob_seed, mymask);

end

