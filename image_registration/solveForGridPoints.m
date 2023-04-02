function mycoordinates = solveForGridPoints(gridPoints)
% it is assumed that our grid points have already been properly aligned! 
%gridPoints should have 6 rows, one for each line. It should have also 4 columns: such
%that the two points establishing the line, (x1,y1, and x2,y2 are set as:
% [x1, y1, x2,y2] 

slopes = zeros(6,1);
yIntercepts = zeros(6,1);

for k = 1:6
    myline = gridPoints(k,:);
slopes(k,1) = (myline(4)-myline(2))/(myline(3)-myline(1)); %rise over run!
yIntercepts(k,1) = myline(2)-slopes(k,1)*myline(1);        % arbitrarily chose my query point to be (x1, y1). plugging that into mx+b=y using the above m to finds b
end


Intersection_points = zeros(2,9);
counter = 1;
for p = 1:3
A = [-slopes(p), 1;-slopes(4), 1];
B = [yIntercepts(p); yIntercepts(4)];
Intersection_points(1:2,counter) = linsolve(A, B);

C = [-slopes(p), 1;-slopes(5), 1];
D = [yIntercepts(p); yIntercepts(5)];
Intersection_points(1:2,counter+1) = linsolve(C, D);

E = [-slopes(p), 1;-slopes(6), 1];
F = [yIntercepts(p); yIntercepts(6)];
Intersection_points(1:2,counter+2) = linsolve(E, F);

counter = counter+3;
end

mycoordinates = Intersection_points';