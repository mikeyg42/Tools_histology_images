function [curveFcn, deriv] = fitTheCurve(x, y, varargin)
% [curveFcn, s] = fitTheCurve(my_img, x, y); 
% OR  [curveFcn, s] = fitTheCurve(my_img, x, y, cornerID, "MOVING/FIXED"); 
% ie:
% varargin{1} = cornerID, the identity of the current corner beging processed
% varargin{2} = identifies wether this is a 'Moving' or a "Fixed'.
% This function tries to fit 5 different non-Rigid transformations o fthe
% control points we had previously established.

% Implementation of the sine curve fitting portion of this script was inspired by the response of
% Star Strider in March 2014 to post on the Mathwworks Q&A page:
% https://www.mathworks.com/matlabcentral/answers/121579-curve-fitting-to-a-sinusoidal-function

% In addition to the sine curve, I also test polynomial curves of order 1, 2
% and 3. Each curve is fit to the data, after which I select the curve maximizing
% R^2 value, indicating that that curve explains a greater proportion of the
% variance of Y. I've also added code for a 4-degree polynomial but one where Ive
% pre-set 2 zero-crossings to be the end points. There might be a mistake here though?...

% Instead of overlaying a grid upon all the data to impose Cartesian space, I
% have opted for a more flexible approach to simplify the curve equations
% being fit, where I position the x-axis to be on the line that intersects the two
% corner points book-ending the domain of each side of a tissue. To achieve
% this, I calculate the angle the aforementioend lines makes with the
% y-axis relative to the entire image, and then rotate the points together
% as a point cloud with an affine transformation about one of the two
% corner points. 

% Michael Glendinning, 2023
warning('off')

%get rid of any contaminating NaN's !!
idx_Not_nan = ~(isnan(x') | isnan(y'));
x = x(idx_Not_nan);
y = y(idx_Not_nan);

% we don't want to use more than 1000 points when fitting our curve!
if numel(x) > 1000
    indx = round(linspace(1, numel(x), 1000)); %linspace will be your way to INDEX the new space
    x = x(indx);
    y = y(indx);
end

yRange = range(y);  
yZeroCrossing = y-mean(y); % optimization is facilitated by the placement of the x axis in the center of the range. 
zeroCrossingIndices = find(diff(sign(yZeroCrossing))~=0);

if any(diff(zeroCrossingIndices)<25)
    idx = find(diff(zeroCrossingIndices)<25);
    zeroCrossingIndices(idx) = round((zeroCrossingIndices(idx)+ zeroCrossingIndices(idx+1))/2);
    zeroCrossingIndices(idx+1) = [];
end

zeroCrossings_x = x(zeroCrossingIndices);
zeroCrossings_x = [min(x); zeroCrossings_x ; max(x)]; % both corners are set as zeroCrossings
zeroCrossings_x = unique(zeroCrossings_x);

% I extract all the combos of point-to-point distances between 0crossings, to pick a suitable average for estimating the period
distances = pdist(zeroCrossings_x);
period = 2*mean(distances(distances>0));% the >0 is necessary because triu leaves a diagonal of zeros

beta0 = [yRange;  period;  0.5;  mean(y)]; % These are the initial "guesses" which the least-squares makes adjutments
sineFit = @(b,x)  b(1).*(sinpi(2*x./b(2) + 2/b(3))) + b(4);    % Format of sine wave fcn to fit
leastSquareCostFcn = @(b) sum((sineFit(b,x) - y).^2);   % This is the Least-Squares cost function to minimize x
[s, ~, ~] = fminsearch(leastSquareCostFcn, beta0); %calling the optimization function FMINSEARCH


% %% fit the polynomial function of 4th degree with the largest and the smallest roots set
% to be the corner points
% --------------------------------------------------
% lets say I have two known roots L and R, and 2 unknown roots a and b:
% y = k*(x-L)(x-R)(x-a)(x-b).... simplifies to:
% y = coeff(1)*
L = min(x(1), x(end));
R = max(x(1), x(end));

polynomFit4 = @(x, coeffs) coeff(1)*(x-L).*(x-R).*(x-coeff(2)).*(x-coeff(3));
ordLeastSquaresCostFcn = @(coeff) sum((polynomFit4(x, coeff) - y).^2);

% add constraints to set L to be the lowest root and R to be the largest root
nonlinear_constraints = @(coeff) deal([], [coeff(2) - L; R - coeff(2); coeff(3) - L; R - coeff(3)]);

% make the initial guesses evenly distributed between L and R
initial_guess = [1, (L+R)/3, (L+R)/1.5]; 

% optimize using fmincon and the SQP algorithm

%opts = optimoptions('fmincon','Algorithm','sqp','Display','none');
opts = optimset(@fmincon, 'Display','none','TolX',1e-6,'TolFun',1e-6);

coeffss = fmincon(ordLeastSquaresCostFcn, initial_guess, [], [], [], [], [], [], nonlinear_constraints, opts);




% --------------------------------------------------
polynomFcn1 = polyfit(x,y, 1);
polynomFcn2 = polyfit(x,y, 2);
polynomFcn3 = polyfit(x,y, 3);

%% my functions
myFncs = {@(x) polyval(polynomFcn1,x), @(x) polyval(polynomFcn2, x), @(x) polyval(polynomFcn3, x), @(x) polyval(p_coeffs_4, x), @(x) sineFit(s,x)};


%% Evaluate 4 curves by calculating the R^2 and RMSE value of each. 
residuals = {y-polyval(polynomFcn1, x), y-polyval(polynomFcn2, x),...
    y-polyval(polynomFcn3, x), y- polyval(p_coeffs_4, x), y-sineFit(s,x)};
resid2 = cellfun(@(x) x.^2, residuals, 'UniformOutput', false);
ssr = cellfun(@(x) sum(x, 'all'), resid2, 'UniformOutput', false);
sst = sum((y-repmat(mean(y), length(y), 1)).^2, 'all');
r_squared = cellfun(@(x) (1-(x/sst)), ssr,  'UniformOutput', false);

rmse = cellfun(@(x) sqrt(mean(x)), resid2,  'UniformOutput', false);

[~, funcIndex_r2] = max(cell2mat(r_squared));
[~, funcIndex_rmse] = min(cell2mat(rmse));

if funcIndex_rmse~=funcIndex_r2
    disp(' inconclusive best fitting function!!');
    funcIndex = max(funcIndex_rmse, funcIndex_r2); 
else
    funcIndex = funcIndex_r2;
end

curveFcn = myFncs{funcIndex};
% we will want the function derivative later on!!! 

if funcIndex<5
    coeffs = {polynomFcn1, polynomFcn1, polynomFcn1,p_coeffs_4};
    deriv = polyder(coeffs{funcIndex});
else
    deriv = @(x) 2*pi*s(1)/s(2).*cospi(2*x./s(2) + 2/s(3)); %this is the derivative with respect to x of our sine curve
end


if isnumeric(varargin{1}) && isstring(varargin{2})
    cornerID = varargin{1};
    imageType = varargin{2};
    if funcIndex < 4
        disp(strcat(num2str(funcIndex), '= degree of polynomial ~~~~~ ',imageType, ' Image, side = ',num2str(cornerID)));
    elseif funcIndex == 4
         disp(strcat('constrained 4th degree polynomial  was used to fit this side ~~~~~ ' , imageType, ' Image, side = ',num2str(cornerID)));
    else
        disp(strcat('sine wave was used to fit this side ~~~~~ ' , imageType, ' Image, side = ',num2str(cornerID)));
    end
else
    disp('issue printing');
end
warning('on')
end
