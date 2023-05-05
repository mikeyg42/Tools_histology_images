
function TF = checkImfill(preFill_bwEdgeMask)
%TF = checkImfill(preFill_bwEdgeMask)
% TF = true then the image passses, it is a good, discontinuity-free edge map.
% TF = false, then it sucks and imfill has ruined everything! 
%
% I usually run this code in the following way:
    % if checkImfill(bwimage)
    %     filledIm = imfill(bwimage)
    % else
    %     something else
    % end
% Michael Glendinning, May 2023

imFilled = imfill(preFill_bwEdgeMask, 'holes');
TF = true;

if sum(imFilled(:)) == sum(preFill_bwEdgeMask(:)) && sum(imFilled(:)) >0 %this means that there were no holes, wit was just a solid mass
disp('imfill had no effect because the binary mask was without any holes already');
return
% this is not a reason to not run imfill.... imfill just won't do anything
end

if dice(bwperim(imFilled), preFill_bwEdgeMask) < 0.5 % @DICE us the Sorensen-Dice similarity measure. 
    % the edge mask and bwperim should correlate >0.99 so if its less than 0.5, no buenoooo
    TF = false;
    return;
end

if sum(imFilled(:)) > 0.93*prod(size(preFill_bwEdgeMask, 1:2)) % undeniable evidence of a leak if this happens
    TF = false;
    return;

elseif sum(imFilled(:)) < 10 || sum(preFill_bwEdgeMask(:)) < 10 % this shouldn't happen...
    disp('input to the @checkImFill function was basically blank !');
    TF = false;
    return;

elseif sum(imFilled(:))/prod(size(preFill_bwEdgeMask, 1:2)) < 0.15 && sum(imFilled(:)) > sum(preFill_bwEdgeMask(:))
    % if sum(imFilled(:)) is not zero but it is less than 15%, that means there were probably
    % some unintended, insignificant tiny skeleton loops filled by 1st call of imfill,

    imFilled = imfill(bwareafilt(imFilled, 1), 'holes');
    if sum(imFilled(:)) > 0.93*prod(size(preFill_bwEdgeMask, 1:2)) || sum(imFilled(:))/prod(size(preFill_bwEdgeMask, 1:2)) < 0.2
        disp('holes inside edge were getting filled. scant changes as a result. try to call bwareafilt(x, 1)...');
        TF = false;
        return;
    end
end
end