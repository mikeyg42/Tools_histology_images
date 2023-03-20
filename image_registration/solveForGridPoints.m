function mycoordinates = solveForGridPoints(gridPoints)

slopes = zeros(6,1);
yIntercepts = zeros(6,1);

for k = 1:6
    myline = gridPoints(k,:);
slopes(k,1) = (myline(4)-myline(2))/(myline(3)-myline(1));
yIntercepts(k,1) = myline(2)-slopes(k,1)*myline(1);
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