function scalingFactor_K = estimateScaling(fixedMat, movingMat)
%fixedMat and movingMat are 4x2 matricies where each row is a different
%corner point. topLeft corner is on top, then it goes clockwise around. 

diag1M =  norm(movingMat(3,:) - movingMat(1,:));
diag2M = norm(movingMat(4,:) - movingMat(2,:));

diag1F =  norm(fixedMat(3,:) - fixedMat(1,:));
diag2F = norm(fixedMat(4,:) - fixedMat(2,:));

% scaling factor "k" is estimated by:
scalingFactor_K = ((diag1F/diag1M) + (diag2F/diag2M))/2; 
end

