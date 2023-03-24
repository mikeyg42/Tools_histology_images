%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  Contents of WSI_processing repo 
% 
%  Michael Glendinning
%  mglendin1@gmail.com
%
%  version 1.0 ~ March 24, 2023
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%  ~~~ IMAGE_SEGMENTATION ~~~
%  backgroundCorrectRGB.m        Background correction (+flatfield correct) using the tophat filter
%  binarizeTissueMG.m            Segmentation of RGB w/ generous threshold on luminosity channel, refined with snakes.      
%  cleanMask.m                   Refine a BW mask by removing blobs too small or too far away from the largest blob
%  ensureDoubleScaled.m          Make sure any amount of images are double floats w/ intensities in range [0,1]
%  fillInPolyLine.m              Redraw a polyline with lots of evenly spaced nodes
%  finalizeMASK.m                Calls and parses responses to the segmentation GUI 
%  geoProbSeg_2tone.m            Wrapper for IPT geodesic segmentation algorithm
%  histeqfloat.m                 Peter Kovesi's histogram equilization function that is  optimized for floating precision
%  hysteresisThresold_w/Morph.m  My implementation of hysteresis to segment 1 blob based on img's histogram and @imextendedmax. 
%  imSegmentation_mosaics.m      START HERE! This function calls iterates through a file dir and segments it semi-automatically
%  LinkUpEdgeDiscontinuities.m   Fill in gaps of an edge map binary image resulting in uninterrupted boundary
%  mikeMedianFilter.m            Median filter with a 5x5 kernel and opt to iterate
%  morphologySegment.m           Segment RGB into binary mask using various morphological   filters
%  preprocessRawRGBims           Calls all the preprocessing functions I like to run on new images
%  segmentationGraphCuts.m       Wrapper function for the IPT function grabcuts to segment foreground
%  textureFilterGUI.m            Opens up a GUI containing many different segmentation options all using textures. Click the best to merge. 
%  vertextTweaking_handdrawn.m   Converts a BW binary image into a polygon, allowing you  to reposition any vertex. Then converts back to BW.
%  weightingPixelIntensityDiffs_RGBseg.m     Segment RGB image by weighting differences in color.

%  ~~~ IMAGE_REGISTRATION ~~~
%  solveForGridPoints.m
%  selectCornersGUI.m            Opens up the GUI to select the 4 corner points 
%  RGB_reorientationGUI.m        Opens up the GUI to coarsely adjust orientation of images at 90degree increments
%  registerSerialSections_part2_nonrigid.m   PART 2 of workflow.
%  registerSerialSections_part1.m            PART 1 of workflow. START HERE!
%  preReg_padding_center.m       Align the centroids of 2 images while also adjusting matricies of any query points using array padding
%  pad4images.m                  Pad set of 2 images+2 masks enough that every white pixel in the mask is some specified distance away from edge.
%  fitTheCurve.m                 Tries to fit 5 different functions (different polynomial or a sine wave) and choose the best using R^2 and the RMSE. 
%  execute3Rotation.m            Executes the rotation as determined by estimateRotation.m
%  estimateScaling.m             Estimates the scaling needed for initial, rigid registration.
%  estimateRotation.m            Uses corner points to estimate how much rotation is needed to align 2 images.
%  curveFittingOfTissueBorders.m Wrapper function for my curve fitting function to sample  each of 4 sides of the tissue piece and fit curve to it
%  createDatastoreWithAll4Images.m            Function pulls the 2 images and their two masks into a datastore to be read into workspace later.
%  arcLength.m                   Calculates the length of an arc. 

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                   
%    Michael Glendinning (mglendin1@gmail.com)
%                       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++