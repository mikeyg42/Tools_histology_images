
function replaceROI(CHOICE, rectums,FINmasks, allROIsamples, outOfBoundsMask, imgCurrent, handy)
close all

allROIsamples(FINmasks{CHOICE}==1) = 0;
FINmasks{CHOICE} = zeros(size(allROIsamples));
rectums(CHOICE).points = [];
xxx = -1;
fHandle = figure;
set(fHandle, 'doublebuffer', 'off')

imDisp = imshow(imgCurrent);

ROISIZE = handy.ROISIZE;

while xxx < 0
    
    win = randomCropWindow2d(size(imgCurrent), ROISIZE);
    rows = win.YLimits(1):win.YLimits(2);
    columns = win.XLimits(1):win.XLimits(2);
    newROI = images.roi.Rectangle(imDisp.Parent, 'Position',[columns(1) rows(1) ROISIZE(2) ROISIZE(1)]);
    maskRect = createMask(newROI);
    
    if sum(maskRect(:)) == 0
        delete(newROI);
        continue
    end
    
    MASKintersect = maskRect.*outOfBoundsMask;
    overlappingROIs = maskRect.*allROIsamples;
    newPoints = newROI.Vertices;
    
    if sum(MASKintersect(:)) == sum(maskRect(:))
        if sum(overlappingROIs(:)) <= 0
            
            visualizeROIs(rectums, allROIsamples, imgCurrent, maskRect, handy);
            
            selection = questdlg('Soooooo?', '~* ROI inquiry *~',...
                'ROI is good','ROI is bad, remake it please', 'ROI is bad, remake it please');
            
            switch selection
                case 'ROI is good'
                    allROIsamples = double(allROIsamples)+double(maskRect);
                    FINmasks{CHOICE} = FINmasks{CHOICE} + maskRect;
                    rectums(CHOICE).points = newPoints;
                    
                    delete(newROI)
                    xxx=47;
                    
                case 'ROI is bad, remake it please'
                    delete(newROI)
                    close all
                    figure;
                    imDisp= imshow(imgCurrent);
                    continue
            end %end of the switch-case block (re: second prompt if new ROI is acceptable
        end %if #3
    end %if #2
    
end % end of the while loop that is responsible for making new ROI and geneerating prompt

showROIs_andthen_pickOne(rectums, FINmasks, allROIsamples, outOfBoundsMask, imgCurrent, handy);
end
