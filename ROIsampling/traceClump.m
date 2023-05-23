function boundary_rc = traceClump(blob_mask)
% boundary = traceClump(blob_mask)
% This function is a wrapper of the builtin function BWTRACEBOUNDARY, automating the
% paramteter selection. Note, this function will only output the EXTERIOR-most boundary.
% Blob=white=foreground

% Michael Glendinning, Apr 2023
% ========================================================================

% find the leftmost point of the blob
leftmost_col = find(sum(blob_mask, 1) > 0, 1, 'first');

% find the highest and lowest points of the blob
[row_inds, ~] = find(blob_mask);
highest_row = min(row_inds);
lowest_row = max(row_inds);

% find the midpoint between the highest and lowest points
mid_row = round((highest_row + lowest_row) / 2);

% get the x, y coordinates of the point we just found
start_col = leftmost_col;
start_row = mid_row;

% use bwperim to get the boundary of the blob
blob_boundary = bwperim(blob_mask, 4);

% use pdist2 to find the point on the boundary closest to our imagined starting point
start_point = [start_row, start_col];
boundary_points = find(blob_boundary(:));
[boundary_rows, boundary_cols] = ind2sub(size(blob_boundary), boundary_points);
boundary_coords = [boundary_rows, boundary_cols];
distances = pdist2(start_point, boundary_coords);
[~, min_idx] = min(distances);
closest_point = boundary_coords(min_idx,:);  % row/col

% Now that we have selected our starting point, we need to determine the starting
% direction (e.g. N, NE, NW, ect.)
fInx = mod(min_idx+3, size(boundary_coords, 1));
if fInx == 0
    fInx =  size(boundary_coords, 1);
end

rInx = mod(min_idx-3, size(boundary_coords, 1));
if rInx == 0
    rInx =  size(boundary_coords, 1);
end

pointFwd = boundary_coords(fInx, :);
pointRev = boundary_coords(rInx, :);
u(:,:) = fliplr(pointFwd-closest_point);
v(:,:) = fliplr(pointRev-closest_point);
zer = [-1, 0];

%Check if u or v is NORTH of closest_point. Use the NORTHERN point to check
if u(:,2)<v(:, 2)
    testangle = abs(atan2(det([u;zer]), dot(u,zer))*180/pi);
else
    testangle = abs(atan2(det([v;zer]), dot(zer,v))*180/pi);
end

%use the angle that northernly point makes with the hoirzonal 
if testangle >135 
  % use bwtraceboundary to trace the boundary of the blob
    boundary_rc = bwtraceboundary(blob_mask, closest_point, 'NE', 8, Inf, 'Clockwise');
elseif testangle < 45 
    boundary_rc = bwtraceboundary(blob_mask, closest_point, 'NW', 8, Inf, 'Clockwise');
else
    boundary_rc = bwtraceboundary(blob_mask, closest_point, 'N', 8, Inf, 'Clockwise');
end

end

