function arcLength = arcLength(deriv, x1, xEnd)
%using the integral of sqrt(1+(dy/dx)^2) formula for arc length

% Define the integrand function
integrand = @(x) sqrt(1 + deriv(x).^2);

% Use quadgk to numerically integrate the integrand over the interval
arcLength = integral(integrand, min(x1, xEnd), max(x1, xEnd));


