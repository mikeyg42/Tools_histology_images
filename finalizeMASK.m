function binaryMask = finalizeMASK(binaryMask,imAdjRGB)
% This function is called refines an initial attempt making mask (binaryMask),
% outlining a particular image (imAdjRGB). Will show overlay of mask on image
% and a GUI with various options to refine the mask, or redraw it entirely. 
% Select DONE to exit

close all force


% we don't need to reredefine this variable each time we loop
listMethods = {'DONE',...
    'Level Set Method segmentation',...
    'Color segment w/ graydiffweight',...
    'Double threshold (hysteresis)',...
    'Color segment with geodesics, w/ provided foreground sample',...
    'Segment w/ Morphology',...
    'k means color clustering segmentation',...
    'Texture segmentation options GUI (no Gabor)',...
    'thresholding / imclose, refined as polyshape',...
    'refine : dilate mask',...
    'refine : just active contour (ie snakes)',...
    'refine : fill in holes',...
    'refine : select only the largest connected blob',...
    'refine : morphology; opening-recon. and closing-recon', ...
    'refine : adjust the vertices of the largest blob',...
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
fig1 = figure; 
axy = axes(fig1);

%initialize visualization
imshow(imAdjRGB,'Parent',axy, 'Border', 'tight');
hold on
visboundaries(axy, binaryMask, 'Color', 'r');
hold off

% The creation of a GUI is within a while loop -- the exit condition for
% which is "MyCounter" becoming positive. This change can onlynhappen when 
% DONE is picked off the menu, causing change in an appdata associated var
% attached to the parent of the figure, ie "the root", the handle of which is 0.

% Until DONE is picked and the appdata feeds a positive value to MyCounter, 
% after a selection is made and cooresponding function runs, the loop
% gets reinitialized. 

drawnow;

MyCounter = -10;
while MyCounter < 0
    
    % this creates the modal dialogue GUI with the list of segementation options.
    [indx, tf] = listdlg('SelectionMode', 'single',...
        'ListString', listMethods,...
        'ListSize', [350, 200]);
    
    if tf == 1
        ImProcessing_choice = listMethods{indx};
        
        if strcmp(ImProcessing_choice,'DONE')
            MyCounter = 10;
            close all force
            
        else
            outOfBoundsMask_old = binaryMask;
            
            close all force
            
            binaryMask = imProcessingResponse(outOfBoundsMask_old, ImProcessing_choice, imAdjRGB);
            
            fig2 = figure('Visible', 'off');
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
        close all force
        disp('GUI closed without selection made. Reopening  now...');
        fig3 = figure;  axy3 = axes(fig3);
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

function outOfBoundsMask = imProcessingResponse(outOfBoundsMask, choice, imgCurrent)

switch choice
   case 'Level Set Method segmentation' 
        outOfBoundsMask = runLevelSet(imgCurrent);
        
    case 'thresholding / imclose, refined as polyshape'
        outOfBoundsMask = binarizeTissueMG(imgCurrent) ;
        
    case 'Double threshold (hysteresis)'
         H = fspecial('gaussian',[3,3],3);
        imgMedian = mikeMedianFilter(imgCurrent, 2, 750, 'RGB');
        current_blurred = imfilter(imgMedian, H);
        imageB = imcomplement(current_blurred);
        grayIm = rgb2gray(imageB); 
        
        outOfBoundsMask = hysteresisThreshold_wMorph(grayIm, imageB);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
        
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
        
    case 'Color segment with geodesics, w/ provided foreground sample'
        outOfBoundsMask = geodesicProbabilitySeg_2colors(imgCurrent);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
    case 'k means color clustering segmentation' 
        outOfBoundsMask = k_means_seg_2colors(imgCurrent,true);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
    case 'refine : morphology; opening-recon. and closing-recon'
        se3 = strel('disk',21);
        Ierode = imerode(outOfBoundsMask,se3);
        Ie_reconstruct = imreconstruct(Ierode, outOfBoundsMask);
        Ir_dilate = imdilate(Ie_reconstruct,se3);
        Iopenbyrecon_closingbyrecon = imreconstruct(imcomplement(Ir_dilate),imcomplement(Ie_reconstruct));
        
        Iopenbyrecon_closingbyrecon = imcomplement(Iopenbyrecon_closingbyrecon);
        outOfBoundsMask = imregionalmax(Iopenbyrecon_closingbyrecon);
        
        outOfBoundsMask= bwareafilt(outOfBoundsMask, 5);
        outOfBoundsMask = cleanMask(outOfBoundsMask);
        
    case 'refine : adjust the vertices of the largest blob'
        outOfBoundsMask  = vertexTweaking_handdrawn(outOfBoundsMask,imgCurrent, 130);
        
    case 'refine : fill in holes'
        outOfBoundsMask = imfill(outOfBoundsMask, 'holes');
        
    case 'refine : select only the largest connected blob'
        outOfBoundsMask = bwareafilt(outOfBoundsMask, 1);
        
    case 'Segment w/ Morphology'
        outOfBoundsMask = morphologySegment(imgCurrent);
end

% SANITY CHECK - some segmentation functions I had to include an imresize call
% because of how slow they would otherwise run...
% and so, just in case I missed resizing, this helped me during debugging
outOfBoundsMaskF = imresize(outOfBoundsMask, size(imgCurrent,1:2));

% similarly this is vestigal from debugging, but gives me peace of mind so its staying!
if ndims(outOfBoundsMaskF)>2
    outOfBoundsMaskF = outOfBoundsMaskF(:,:,1);
end

end %go back to the while loop in the above function (@finalizeMask)