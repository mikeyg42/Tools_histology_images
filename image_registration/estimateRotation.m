function [theta] = estimateRotation(movingMat,fixedMat)
%syntax = [theta] = estimateRotation(movingMat,fixedMat)
% This is an ESTIMATE! I calculate the angles formed by the lines 
% connecting the top or bottom two corners of
% the mask with the y axis. I take the mean of my 2 estimates for fixed
% and 2 estimates for the moving, 
% then subtract: fixed - moving = dTheta, which is used to generate rotMatrix

vec_down = [0,-1];

vecM_acrossB = movingMat(3,:)-movingMat(4,:);
vecM_acrossT = movingMat(2,:)-movingMat(1,:);

vecF_acrossB = fixedMat(3,:)-fixedMat(4,:);
vecF_acrossT = fixedMat(2,:)-fixedMat(1,:);

thetaMbot = acosd(dot(vec_down,vecM_acrossB)/(norm(vec_down)*norm(vecM_acrossB)));
thetaMtop = acosd(dot(vec_down,vecM_acrossT)/(norm(vec_down)*norm(vecM_acrossT)));
delta_moving = (thetaMbot+thetaMtop)/2;

thetaFbot = acosd(dot(vec_down,vecF_acrossB)/(norm(vec_down)*norm(vecF_acrossB)));
thetaFtop = acosd(dot(vec_down,vecF_acrossT)/(norm(vec_down)*norm(vecF_acrossT)));
delta_fixed = (thetaFbot+thetaFtop)/2;

theta = delta_fixed-delta_moving; 

