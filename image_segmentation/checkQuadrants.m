
function TF = checkQuadrants(edgeIm)
%TF = checkQuadrants(edgeIm)
% Check if there is but one tiny speck of white in each quadrant of your image.
% if there is 1 quadrant with zero signal, the function will return back a FALSE or 0 for
% BAD. if everything is GOOD and there is some white everywhere, the function returns a
% TRUE or a 1.  1 = good, 0 = bad.
% Michael Glendinning, Apr 2023
% ============================================================================================
% confirm that you have either a logical image OR a float matrix with only 0s and 1s
if ~islogical(edgeIm) && ~(all(unique(edgeIm(:))==0|unique(edgeIm(:))==1) && isfloat(edgeIm))
    edgeIm = logical(edgeIm); %all nonzero values become 1 and zeros become 0
    disp('checkQuadrant function expected either a logical mat or a float with only 0s and 1s');
end
sz = size(edgeIm, 1:2);
sz2 = sz./2;
tf = any(sum(edgeIm( 1:ceil(sz2(1))     ,1:ceil(sz2(2))     ), 'all') == 0 |...
         sum(edgeIm( floor(sz2(1)):sz(1),1:ceil(sz2(2))     ), 'all') == 0 | ...
         sum(edgeIm( 1:ceil(sz2(1))     ,floor(sz2(2)):sz(2)), 'all') == 0 |...
         sum(edgeIm( floor(sz2(1)):sz(1),floor(sz2(2)):sz(2)), 'all') == 0);
TF = ~tf;
end