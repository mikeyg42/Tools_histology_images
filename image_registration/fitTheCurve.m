function [curveFcn, derivCoeffs] = fitTheCurve(x, yP, cornerID, varargin)
% [curveFcn, s] = fitTheCurve(my_img, x, y, cornerID, ...
% {optional: an identifierImage}, {optional: verbose_y/n}); 
% ie:
% cornerID = the identity of the current corner beging processed
% varargin{1} = identifies wether this is a 'Moving' or a "Fixed'. only for the read out
% at the end.
% varargin{2} = verbose readout true or false (or 0 or 1). default is verbose ON

% This function tries to fit a lot of polynomial curves to the data. It also tries a sine
% curve. Please keep in mind that we have no reason to worry about "over-fitting".... in
% fact, we WANT to overfit. This is because we are not trying to predict anything with our
% curve... we just want a very accurate, differentiable function
% easily and it will allow is to place control points on a jagged edge. 

% Implementation of the sine curve fitting portion of this script was inspired by the response of
% Star Strider in March 2014 to post on the Mathwworks Q&A page:
% https://www.mathworks.com/matlabcentral/answers/121579-curve-fitting-to-a-sinusoidal-function
% I've added to it considerably

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
%========================================================================

assert(isnumeric(cornerID)); %should be 1, 2, 3, 4

if isempty(varargin{1})
    imageidentity = 'image'; % a generic filler identifier
else
    imageidentity = varargin{1};
end

try assert(isstring(imageidentity))
catch
    assert(ischar(imageidentity))
end
imageType = imageidentity;

try  verbose_tf = varargin{2};
catch
    verbose_tf = true;
end

%get rid of any contaminating NaN's !!
idx_Not_nan = ~(isnan(x') | isnan(yP'));
x = x(idx_Not_nan);
yP = yP(idx_Not_nan);

% now we will bring the range down to be centered on the y-axis! we will bring it back up
% to original level at the end
mn = mean(yP);
y = yP-mn;

xunsamp = x;
yunsamp = y;
sampSize = {100, 250, 400, 800};
kDeg = [1:1:15, 12:3:24, 4];% degrees of polynomial
s = struct([]);
for k = 1:size(sampSize, 2)
    indx = unique(round(linspace(1, numel(xunsamp), sampSize{k})));
    x = xunsamp(indx);
    y = yunsamp(indx);
    
    s(k).pFcn1 = polyfit(x,y, kDeg(1));
    s(k).pFcn2 = polyfit(x,y, kDeg(2));
    s(k).pFcn3 = polyfit(x,y, kDeg(3));
    s(k).pFcn4 = polyfit(x,y, kDeg(4));
    s(k).pFcn5 = polyfit(x,y, kDeg(5));
    s(k).pFcn6 = polyfit(x,y, kDeg(6));
    s(k).pFcn7 = polyfit(x,y, kDeg(7));
    s(k).pFcn8 = polyfit(x,y, kDeg(8));
    s(k).pFcn9 = polyfit(x,y, kDeg(9));
    s(k).pFcn10 = polyfit(x,y,kDeg(10));
    s(k).pFcn11 = polyfit(x,y,kDeg(11));
    s(k).pFcn12 = polyfit(x,y,kDeg(12));
    s(k).pFcn13 = polyfit(x,y,kDeg(13));
    s(k).pFcn14 = polyfit(x,y,kDeg(14));
    s(k).pFcn15 = polyfit(x,y,kDeg(15));
    
    s(k).ransac12 = fitPolynomialRANSAC([x,y],kDeg(16),20);
    s(k).ransac15 = fitPolynomialRANSAC([x,y],kDeg(17),20);
    s(k).ransac18 = fitPolynomialRANSAC([x,y],kDeg(18),20);
    s(k).ransac21 = fitPolynomialRANSAC([x,y],kDeg(19),20);
    s(k).ransac24 = fitPolynomialRANSAC([x,y],kDeg(20),20);
    
    [sines{k}, cof{k}] = fitSineCurve(x,y);
end

nFuncs = numel(fieldnames(s))+1;
nSampSizes = size(sampSize, 2);

coeffs = struct2cell(s); %extract coefficients from your strucutral aray
yest = cellfun(@(k) polyval(k,x), coeffs, 'UniformOutput', false); %use coeffs as coeffs in functions to evaluate y estimate

yest{nFuncs, 1, 1} = sines{1}(x); %estimate y with the sine curve and incorporate into yEsst cell arrya
yest{nFuncs, 1, 2} = sines{2}(x);
yest{nFuncs, 1, 3} = sines{3}(x);
yest{nFuncs, 1, 4} = sines{4}(x);

resid = cellfun(@(f) y-f, yest, 'UniformOutput', false); %calculate residuals
residSq = cellfun(@(x) x.^2, resid, 'UniformOutput', false);% take the square of the residuals

%calculate R^2 using SSR and SST
ssr = cellfun(@(x) sum(x, 'all'), residSq, 'UniformOutput', false);
sst = sum((y-repmat(mean(y), length(y), 1)).^2, 'all');
rsq = cellfun(@(x) (1-(x/sst)), ssr,  'UniformOutput', false);

%this is R^2 
rsq = squeeze(cell2mat(rsq)); 

%calculate adjusted R^2
[Sidx,~] = meshgrid(1:size(sampSize,2), 1:nFuncs);
nSamps = reshape([sampSize{Sidx}], [], nSampSizes);
kParams = repmat(kDeg',1, nSampSizes);
oo = ones(nFuncs,nSampSizes);
adjR_squared = oo-((oo - rsq).*(nSamps - oo)./(nSamps - kParams - oo));

%find the highest R^2 value
[rsquared_win, fIx] = max(adjR_squared, [], 'all', 'linear');

%% rewrite function into an anonymous function format 
%(... would be much easier with symbolic math toolbox function poly2sym)

bestRow = mod(fIx, nFuncs);
samp = ceil(fIx/nFuncs); %indicates which sample size made best function
if bestRow ~= nFuncs %this will be true only for when the best fitting curve is the sine curve function
    
    allfnames = fieldnames(s); % to retrieve the function coeffs we need to index using the field name
    fName = allfnames(bestRow);
    coeffs = getfield(s, {samp}, char(fName));
    fullcoeffs = zeros(1, 1+max(kDeg)-size(coeffs, 2));
    curveFcn = [fullcoeffs, coeffs];
    curveFcn(1, max(kDeg)+1) = curveFcn(1, max(kDeg)+1)+mn; %add back in mean(y) that we subtracted away
    %curveFcn = @(x) ([x.^24, x.^23, x.^22, x.^21, x.^20, x.^19, x.^18, x.^17, x.^16, x.^15, x.^14, x.^13, x.^12, x.^11, x.^10, x.^9, x.^8, x.^7, x.^6, x.^5, x.^4, x.^3, x.^2, x.^1, x.^0].*coeffsFull);
    
    %% calculate the derivative:
    d = kDeg(bestRow); %degree of your polynomial
    
    derivCoeffs = coeffs(1:end-1).*[d:-1:1]; % This is how the coefficients are calculated...
    
    %fullD = zeros(1, max(kDeg)-size(derivCoeffs, 2));
    %derivFull = [fullD, derivCoeffs];
    %myDeriv = @(x) ([x.^23, x.^22, x.^21, x.^20, x.^19, x.^18, x.^17, x.^16, x.^15, x.^14, x.^13, x.^12, x.^11, x.^10, x.^9, x.^8, x.^7, x.^6, x.^5, x.^4, x.^3, x.^2, x.^1, x.^0].*derivFull);
    
else %this is just for the sine curve
    g = cof{samp};
    g(4) = g(4)+mn; %add back in the mean(y) we subtracted away earlier!
    curveFcn = @(x) g(1).*(sinpi(2*x./g(2) + 2/g(3))) + g(4)+ g(5).*(sinpi(2*x./g(6) + 2/g(7)));
    %myDeriv = @(x) 2*pi*(g(1)/g(2).*cospi(2*x./g(2) + 2/g(3))+ g(5)/g(6).*(cospi(2*x./g(6)+ 2/g(7))));

end
    
% turn warnings back on 
warning('on');

%% verbose readout
if verbose_tf %explicitly indicating a zero here means the user wants no verbosity
    
    if bestRow == nFuncs %ie sine wave
        disp(strcat('sine waveform was best fit 2 this side! Adj R^2=',num2str(rsquared_win),'type of img: ', imageType, ' sideID = ',num2str(cornerID))); 
    
    elseif contains(fName, 'ransac') %gets the RANSAC polys
        disp(strcat('RANSAC polynomial of deg:',num2str(d),'. Adj R^2 = ', num2str(rsquared_win), 'And, type of image = ', imageType, 'sideID = ',num2str(cornerID)));
    
    else
        disp(strcat('Polyfit, deg = ',num2str(d),'. Adj R^2 = ', num2str(rsquared_win), 'And, type of image = ', imageType, 'sideID = ',num2str(cornerID)));
        
    end
end

end


