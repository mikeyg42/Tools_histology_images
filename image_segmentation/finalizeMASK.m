function binaryMask = finalizeMASK(mySettings, binaryMask,imAdjRGB)
% This function is called refines an initial attempt making mask (binaryMask),
% outlining a particular image (imAdjRGB). Will show overlay of mask on image
% and a GUI with various options to refine the mask, or redraw it entirely. 
% Select DONE to exit

%Michael Glendinning, 2023
% ======================================================================== %

close all force

% we don't need to reredefine this variable each time we loop
listMethods = {'DONE',...
    'Simple thresholding, plus some gentle cleaning',...
    'Color segment w/ graydiffweight',...
    'Segment with various vector calc tools'
    'Full color segmentation',...
    'Double threshold (hysteresis)',...
    'Color segment with geodesics, w/ provided foreground sample',...
    'Segment w/ Morphological transformations',...
    %'Level Set Method segmentation',...
    'k means color clustering segmentation',...
    'Texture segmentation options GUI (no Gabor)',...
    %'(needs work still) Texture segmentation Gabor+LBP, then FCM',...
    'refine : dilate mask',...
    'refine : just active contour (ie snakes)',...
    'refine : fill in holes',...
    'refine : select only the largest connected blob',...
    'refine : morphology; opening-recon. and closing-recon',...
    'refine : adjust the vertices of the largest blob',...
    'refine : multiple blobs? check for outliers',...
    'NOT DONE! Save as is and close',...
    };  
 
% Sanity Check - make sure mask and image are the same size!
size1 = size(binaryMask, 1:2);
size2 = size(imAdjRGB, 1:2);
if size1 ~=size2
    disp('unequal sized images! resizing the mask!');
    binaryMask = imresize(binaryMask, size2, {@oscRsesampling,4});
end 
%--------------- Make the list dialogue gui loop -----------------
% save time by making image display figure outside the loop
fig1 = uifigure('CloseRequestFcn',@guiCloseFunction); 
axy = axes(fig1);

%initialize visualization
imshow(imAdjRGB,'Parent',axy, 'Border', 'tight');
hold on
visboundaries(axy, binaryMask, 'Color', 'r');
hold off

% The creation of a GUI is within a while loop -- the exit condition for
% which is "MyCounter" becoming positive. This change can onlynhappen when 
% DONE is picked off the menu. 

% Until DONE is picked and the appdata feeds a positive value to MyCounter, 
% after a selection is made and cooresponding function runs, the loop
% gets reinitialized. 

drawnow;

while false
    
    % this creates the modal dialogue GUI with the list of segementation options.
    [indx, tf] = listdlg('SelectionMode', 'single',...
        'ListString', listMethods,...
        'ListSize', [360, 210],...
        'CancelString', 'Quit without saving');
    
    if tf == 1
        ImProcessing_choice = listMethods{indx};
        close all force
        if strcmp(ImProcessing_choice,'DONE')
            break

        else
            outOfBoundsMask_old = binaryMask;
            
            binaryMask = imProcessingResponse(mySettings, outOfBoundsMask_old, ImProcessing_choice, imAdjRGB);
            
            fig2 = uifigure('Visible', 'off', 'CloseRequestFcn',@guiCloseFuntion);
            axy2 = axes(fig2);
            
            %show new boundary in yellow, old boundary in red
            imshow(imAdjRGB,'Parent',axy2, 'Border', 'tight');
            hold on
            visboundaries(axy2, outOfBoundsMask_old, 'Color', 'r');
            visboundaries(axy2, binaryMask, 'Color', 'b');
            hold off
            
            fig2.Visible = 'on';
            drawnow;
        end
    elseif tf == 0 
% tf=0 when: 
        % - 'x'-out of the GUI or use ESC to close window
        % - 'cancel' button is clicked
        % - 'not done! save and quit' is selected  off the list
        close all force
        
        %recreate the figure
        fig3 = uifigure('CloseRequestFcn',@guiCloseFuntion);  
        axy3 = axes(fig3);
        
        imshow(imAdjRGB,'Parent',axy3, 'Border', 'tight');
        hold on
        visboundaries(axy3, binaryMask, 'Color', 'r'); hold off
    end
end

close all force
end

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  %
% ^    ^      ^     ^   This function creates the list menu GUI within a while-loop 
% |                 |   that keeps reappearing after each fcn call, until DONE is selected. 
% |___The GUI Fcn __|   This also shows a visualization of results of each Fcn.
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  %

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  %
% |---Callbacks----|   This is a "relay" function that connects user selection 
% |   for GUI      |   to corresponding function, each of which has
% |                |      INPUT = imgCurrent (class double;RGB)
% |                |   (+ for "refine" fcns, there is a 2nd input of the draft mask)
% V                V      OUTPUT = outOfBoundsMask (class binary; foregr(1)/backgr(0))
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  %

function outOfBoundsMask = imProcessingResponse(mySettings, outOfBoundsMask, choice, imgCurrent)

switch choice
        
    case 'Simple thresholding, plus some gentle cleaning'
        outOfBoundsMask = binarizeTissueMG(imgCurrent);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
    case 'Double threshold (hysteresis)'
         H = fspecial('gaussian',[3,3],3);
        imgMedian = mikeMedianFilter(imgCurrent, 2, 750, 'RGB');
        current_blurred = imfilter(imgMedian, H);
        imageB = imcomplement(current_blurred);
        
        outOfBoundsMask = hysteresisThreshold_wMorph(imageB);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
    
    case 'Full color segmentation'
        imgCurrent = imresize(imgCurrent, 0.5, {@oscResampling, 4});
        [BWmask, ~] = fullColorSegmentation(imgCurrent);
        BWmask = imresize(BWmask, 2, {@oscResampling, 4});
    
    case 'Segment with various vector calc tools'
        myIm_resized = imresize(imgCurrent, round(5000/max(size(imgCurrent, 1:2)), 3) , {@oscResampling, 4} );
        outOfBoundsMask = vectorCalcSegment(myIm_resized);
        outOfBoundsMask = imresize(outOfBoundsMask, size(imgCurrent, 1:2) , {@oscResampling, 4});

    case 'Color segment w/ graydiffweight' 
        
        outOfBoundsMask = weightingPixelIntensityDiffs_RGBseg(imgCurrent);
         
    case 'refine : dilate mask'
        outOfBoundsMask2 = bwpack(outOfBoundsMask);
        se3 = strel('disk', 20, 8);
        outOfBoundsMask2 = imdilate(outOfBoundsMask2, se3, 'ispacked');
        outOfBoundsMask = bwunpack(outOfBoundsMask2, size(outOfBoundsMask,1));
        
        imgCurrent_gray = rgb2gray(imgCurrent);
        outOfBoundsMask = activecontour(imgCurrent_gray, outOfBoundsMask, 60, 'edge', 'ContractionBias', 0.9);
        
    case 'refine : just active contour (ie snakes)'
        imgCurrent_gray = rgb2gray(imgCurrent);
        imgCurrent_gray = imresize(imgCurrent_gray, 0.5, {@oscResampling, 4});
        imgFin = mikeMedianFilter(imgCurrent_gray, 2, 1250, 'Grayscale'); %2 iterations/ downsampled so largest side becomes 1250. reverts before snakes
        outOfBoundsMask = activecontour(imgFin, outOfBoundsMask, 75, 'Chan-Vese', 'SmoothFactor' ,0.15);
        outOfBoundsMask = imresize(outOfBoundsMask, 2, {@oscResampling, 4});
        
    case 'Texture segmentation options GUI (no Gabor)' 
        rescaleFactor = 2000/max(size(imgCurrent, 1:2));
        imgRGB = imresize(imgCurrent,rescaleFactor, {@oscResampling,4}); %gets resized after the switch-case block 
        imgGray = rgb2gray(imgRGB);
        outOfBoundsMask = textureFilterGUI(imgGray); 
        outOfBoundsMask = imresize(outOfBoundsMask, size(imgCurrent, 1:2),{@oscResampling,4});
        
    case 'Color segment with geodesics, w/ provided foreground sample'
        outOfBoundsMask = geoProbSeg_2tone(imgCurrent);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
    case 'k means color clustering segmentation' 
        outOfBoundsMask = k_means_seg_2colors(imgCurrent,true);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
    case '(needs work still) Texture segmentation Gabor+LBP, then FCM'
        gray_img = rgb2gray(imgCurrent);
        outOfBoundsMask = textureSeg_FCM(gray_img);
    
    case 'Level Set Method segmentation' 
        outOfBoundsMask = runLevelSet(imgCurrent);

    case 'refine : morphology; opening-recon. and closing-recon'
        se3 = strel('disk',21);
        Ierode = imerode(outOfBoundsMask,se3);
        Ie_reconstruct = imreconstruct(Ierode, outOfBoundsMask);
        Ir_dilate = imdilate(Ie_reconstruct,se3);
        Iopenbyrecon_closingbyrecon = imreconstruct(imcomplement(Ir_dilate),imcomplement(Ie_reconstruct));
        
        Iopenbyrecon_closingbyrecon = imcomplement(Iopenbyrecon_closingbyrecon);
        outOfBoundsMask = imregionalmax(Iopenbyrecon_closingbyrecon);
        
    case 'refine : multiple blobs? check for outliers'
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
    case 'refine : adjust the vertices of the largest blob'
        outOfBoundsMask  = vertexTweaking_handdrawn(outOfBoundsMask,imgCurrent, 130);
        
    case 'refine : fill in holes'
        outOfBoundsMask = imfill(outOfBoundsMask, 'holes');
        
    case 'refine : select only the largest connected blob'
        outOfBoundsMask = bwareafilt(outOfBoundsMask, 1);
        
    case 'Segment w/ Morphological transformations'
        outOfBoundsMask = morphologySegment(imgCurrent);
    
    case 'NOT DONE! Save as is and close'
        [filepath, ~, ~] = get_filepath(mySettings.directories.rawData, mySettings.activeFilename);
        info = imfinfo({filepath});
        segment_SaveFunc(mySettings, info, outOfBoundsMask, imgCurrent);
        close all force
        drawnow 

        error('Not an error, just an update: Execution stopped by the user. Files were saved and progress can be resumed at a later time. Cheers!');

end

% sanity check #1
outOfBoundsMaskF = imresize(outOfBoundsMask, size(imgCurrent,1:2));

%sanity check #2
if ndims(outOfBoundsMaskF)>2
    outOfBoundsMaskF = outOfBoundsMaskF(:,:,1);
end

end %end improcessing function. returns now to parent function workspace (@finalizeMask))



function guiCloseFuntion(src, ~, varargin)
% guiCloseFuntion(src, evt, varargin)
% 
% This is a callback function for the GUIs to prevent bugs associated with inadvertent
% closing of GUIs and/or mis-selection of a particular choice. 
if nargin>2
    fHandle = varargin{1};
else
    fHandle = src;
end

mymessage = sprintf(['<html> <font size="5"]><b>Are your sure you would like to close right now? without saving?!</b></font<br><br>' ...
'<font size="3">You can save by selecting <font color="blue"><i>DONE</i></font> or <font color="blue"><i>NOT DONE</i></font> on the list to initiate this. </font></html>']);  

confirmation = uiconfirm(fHandle, mymessage,'Confirm Close',...
    'Options', {'Yes, get me out of here, no saving', 'No, Ive changed my mind'},...
    'DefaultOption', 'No, Ive changed my mind',...
    'Icon', 'warning' ,...
    'CancelOption', 'No, Ive changed my mind',...
    'Interpreter', 'html');

switch confirmation
    case 'Yes, get me out of here, no saving'
       
        close all force
        error('Not really an error, just an update: Execution stopped by the user and, as requested, not saved. Cheers!');

    case 'No, Ive changed my mind'
        return;
end
end

 
