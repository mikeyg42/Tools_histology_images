function [sineFunc1, h] = fitSineCurve(x,y)
% [sineFunc1, coeffs] = fitSineCurve(x,y)

options = optimset('Display','off', 'MaxIter',1500);

yZeroCrossing = y-mean(y); % optimization is facilitated by the placement of the x axis in the center of the range.
zeroCrossingIndices = find(diff(sign(yZeroCrossing))~=0);

if any(diff(zeroCrossingIndices)<15) %this confirms that my zero crossings aren't too close to each other
    idx = find(diff(zeroCrossingIndices)<15);
    zeroCrossingIndices(idx) = round((zeroCrossingIndices(idx)+ zeroCrossingIndices(idx+1))/2);
    zeroCrossingIndices(idx+1) = [];
end

zeroCrossings_x = x(zeroCrossingIndices);
zeroCrossings_x = [min(x); zeroCrossings_x ; max(x)];
zeroCrossings_x = unique(zeroCrossings_x);

% I extract all the combos of point-to-point distances between 0crossings, to pick a suitable average for estimating the period
distancesx = pdist(zeroCrossings_x);
periodx = 2*mean(distancesx(distancesx>0));% the >0 is necessary because triu leaves a diagonal of zeros

beta0 = [range(y);  periodx;  0.5;  mean(y)]; % These are the initial "guesses" which the least-squares makes adjutments
sineFit = @(b,x)  b(1).*(sinpi(2*x./b(2) + 2/b(3))) + b(4);    % Format of sine wave fcn to fit
leastSquareCostFcn = @(b) sum((sineFit(b,x) - y).^2);   % This is the Least-Squares cost function to minimize x
[s, ~, exitflag] = fminsearch(leastSquareCostFcn, beta0, options); %calling the optimization function FMINSEARCH
sineFunc1 = @(x) s(1).*(sinpi(2*x./s(2) + 2/s(3))) + s(4);
%%
if exitflag == 0
    beta1 = (beta0+s)./2;
    [s2, ~, ~] = fminsearch(leastSquareCostFcn, beta1, options);
    sineFunc2 = @(x) s2(1).*(sinpi(2*x./s2(2) + 2/s2(3))) + s2(4);
    
    res = {y-sineFunc1(x), y-sineFunc2(x)};
    resid2 = cellfun(@(x) x.^2, res, 'UniformOutput', false);
    ssr = cellfun(@(x) sum(x, 'all'), resid2, 'UniformOutput', false);
    sst = sum((y-repmat(mean(y), length(y), 1)).^2, 'all');
    r_squared = cell2mat(cellfun(@(x) (1-(x/sst)), ssr,  'UniformOutput', false));
    if r_squared(1)<r_squared(2)
        sineFunc1 = sineFunc2;
        s = s2;
    end
end

sineFit2 = @(c, x) sineFunc1(x) + c(1).*(sinpi(2*x./c(2) + 2/c(3))); %I changed the sign of the term within sine to make waves more out of phase
beta0beta = [6;  3;  0.5];
leastSquareCostFcn3 = @(c) sum((sineFit2(c,x) - y).^2);
[w, ~, ~] = fminsearch(leastSquareCostFcn3, beta0beta, options); %calling the optimization function FMINSEARCH
h = [s(1), s(2), s(3), s(4), w(1), w(2), w(3)];

sineFuncFull = @(x) h(1).*(sinpi(2*x./h(2) + 2/h(3))) + h(4)+ h(5).*(sinpi(2*x./h(6) + 2/h(7)));

res = {y-sineFunc1(x), y-sineFuncFull(x)};
resid2 = cellfun(@(x) x.^2, res, 'UniformOutput', false);
ssr = cellfun(@(x) sum(x, 'all'), resid2, 'UniformOutput', false);
sst = sum((y-repmat(mean(y), length(y), 1)).^2, 'all');
r_squared = cell2mat(cellfun(@(x) (1-(x/sst)), ssr,  'UniformOutput', false));
if r_squared(1)<r_squared(2)
    sineFunc1 = sineFuncFull;
    h(5) = 0;
end


end

