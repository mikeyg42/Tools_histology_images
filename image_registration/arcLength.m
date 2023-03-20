function arcLen = arcLength(xQueryPoints, yValues)
%using the sqrt(1+(dy/dx)^2) formula for arc length but dy/dx is calculated

dydxv=diff(yValues)./diff(xQueryPoints); %linear spline between each node (this vectorized form of y(i+1)-yi divided by x(i+1)-xi

m = length(xQueryPoints);
dydxv(m) = dydxv(m-1);

integrand = sqrt(1+dydxv.*dydxv);
arcLen = trapz(xQueryPoints,integrand);
    end
