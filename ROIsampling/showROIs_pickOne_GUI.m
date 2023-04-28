
function showROIs_pickOne_GUI(ROIcentroidCoords, FINmasks, BWimageOfROIS, outOfBoundsMask, imgCurrent, handy)
%% next we ask if we should exit and move to another image, or replace one?
NUMROI = handy.NUMROI;

visualizeROIs(ROIcentroidCoords, BWimageOfROIS, imgCurrent, [], handy)

numlist = num2cell(1:NUMROI);
exitChoice = {'all good!'};

CHOICE = menu('Toss out of these any ROIs?', [numlist, exitChoice]);
exitNumber = NUMROI+1;

switch CHOICE
    case exitNumber
        
        disp(strcat('DONE with:  ', handy.filenames{handy.counter}));
        
        saveFcn(FINmasks, BWimageOfROIS, outOfBoundsMask, handy)
        
    case numlist
        
        replaceROI(CHOICE, ROIcentroidCoords,FINmasks, BWimageOfROIS, outOfBoundsMask, imgCurrent, handy);
        
end
end
