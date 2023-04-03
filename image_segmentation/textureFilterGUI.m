function finMask = textureFilterGUI(imgC)
%% This function is a reliable semi-automated segmentation tool. 
% Grayscale images are entered, and then 8 different possible binary
% masks are generated. Most of these segmentation techniques are classic
% in the field and implemented without substantial conceptual modifcation. In various ways
% I set parameters optimally for my data, so if poor performance in encountered, you will
% want to run through the code for those instances. 
% There are 8 images to pick from:
%   - 1. image gradient, and is implemented based on Matlab documentation
%   - 2. some manipulation of the first image to
%   ensure that I've segemented the largest 2.5% of blobs in the image. 
%   - 3, 4, 5. all rely on different means of locally filtering 
% in the spatial domain so to assess regional differences in specific textures.
% I've kept the convolutions here at the recommended 9x9 square kernel. I
% have entropy(areas exhibiting randomness are brighter), range (brighter =
% more range surrounding a particular pixel), and std of course has
% brighter pixels indicative of more variation in neighborhood
% .  - Next I have the very classical gray-level co-occurance matrix
% derived image
%    - After that there is commmented out the option to filter using 
%    matrix factorization to detect textures. It works well, just was too
%    slow with my large images. I rewrote to hasten, but still wasn't
%    adding much to the party. 
% see : https://www.mathworks.com/matlabcentral/fileexchange/48817-svd-eigen-qr-and-lu-texture-transforms
%    - imextendedmax and min use the distance transform between known
%    boundary points, I think at least, and then converts to H-maxima transforms.
%    similar to the watershed seg. algorigthm
%   
%    

        se1 = strel('diamond',7);
        se2 = strel([0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0]);
        se3 = strel('disk', 13);
    % Nice classical way to increase contrast quickly 
        imgC = imsubtract(imadd(imgC,imtophat(imgC,se1)),imbothat(imgC,se1));

    % gradient image filtering + blob analysis
        [Gx,Gy] = imgradientxy(imgC, 'intermediate');
        [Gmag,~] = imgradient(Gx,Gy);
        BW1 = imclose(imbinarize(rescale(Gmag, 0, 1),'global'), se3);
        BW1 = imfill(BW1, 'holes');% maybe better to do: imcomplement(imfill(imcomplement
        
   %blob analysis for the Gmag image BW1   
        [labxIm, numBlobs] = bwlabel(BW1);
        blobAreas = regionprops(labxIm, 'area');
        allAreas = [blobAreas.Area];
        [~, IdxSort] = sort(allAreas, 'descend');
        biggestBlob = ismember(labxIm, IdxSort(1:ceil(numBlobs*0.025))); %takes top 2.5%
        BW2 = biggestBlob > 0;
        BW2 = imfill(imclose(BW2, se2), 'holes');
        
    % texture filters: entropy, range, variance       
        Eim = rescale(entropyfilt(imgC));
        Sim = rescale(stdfilt(imgC,ones(9)));
        Rim = rescale(rangefilt(imgC,ones(9)));
        imsTexture = {Eim, Sim, Rim};
        imsBinarized = cellfun(@(x) imbinarize(x, 'global'), imsTexture, 'UniformOutput', false);
        imsBinarized = cellfun(@(x) imclose(x, se2), imsBinarized, 'UniformOutput', false);
        imsBinarized = cellfun(@(x) imfill(x,'holes'), imsBinarized, 'UniformOutput', false);
        BW3 = imsBinarized{1};
        BW4 = imsBinarized{2};
        BW5 = imsBinarized{3};
    % gray-level co-occurance matrix - 
% the GLCM checks each pixel and evaluates the extent to which the intensity of a pixel is shared 
% by pixels its spatially adjacent neighbors. See:  
% Haralick, R.M. et al,(1997) "Textural Features for Image Classification" 

% the offset mnatrix can be made bigger or smaller! Larger risks decreased resolution and worsened boarder effects. 

        offsets = [0 1; -1 1;-1 0;-1 -1]; 
        % ^ this basically 8-point connectivity (here says NW, N, NE, E.. 
        [~,ScaledI] = graycomatrix(imgC,'Offset',offsets, 'Symmetric', true);
        % ...and with symmetric indicated true, we include the remaining other 4/8 directions
        BW6 = imbinarize(rescale(imcomplement(ScaledI)), 'global');  
     
     % Lu matrix texture transform - effective but TOO SLOWWWWWWWW! forked
     % from online - maybe I can make it faster in the future?
     
     %   BW7 = textureTransform(imgC); 
     
     % maximum AND minimum regions(not explicitly a texture filter... but
     % they are effective here and conceptually adjacent)
     
     Itop = imtophat(imgC, strel('diamond', 15));
     BW7a = imextendedmin(Itop, 0.5, 8);
     BW7a = bwareafilt(BW7a, 1);
     BW7 = imcomplement(BW7a); %minimum image
     
     BWa = imextendedmax(Itop, 0.1, 8);
     BW8 = imclose(BWa, strel('disk', 9, 8));
     
     % last step quickly draw an arbitrary background mask - just as a sanity check, 
     % I double check that the edge band (ostensibly all background) is not white. 
     % If it were (ie mean of the region>0.5), it is inverted.
     SZ = size(BW8, 1:2);
     mask = false(SZ);
     mask(1:3,1:SZ(2))=true; 
     tester8 = mean2(BW8(mask));
     tester7 = mean2(BW7(mask));
     if tester8 > 0.5, BW8 = imcomplement(BW8); end
     if tester7 > 0.5, BW7 = imcomplement(BW7); end

%% PART 2: images are made = time to select the winners and combine their wisdom
     % 
     ImList = {BW1,BW2,BW3,BW4 BW5, BW6, BW7,BW8};
     out = imtile(ImList, 'GridSize', [2, 4],...
         'BorderSize', 5);
     
     figy = uifigure('Visible', 'off');
     gl_master = uigridlayout(figy, [5, 1]);
     bg = uibuttongroup(gl_master, 'Visible', 'on');
     
     pnl = uipanel(gl_master, 'Visible', 'on');
     gl3_images = uigridlayout(pnl, [7, 1]);
     ax = uiaxes(gl3_images);
     ax.GridLineStyle = 'none';
     imshow(out, 'Parent', ax);
     disableDefaultInteractivity(ax);
     
     pb1= uibutton(gl3_images,'Text', 'Confirm all good masks are checked, so we can move on',...
         'ButtonPushedFcn',@buttonCallback);

     gl2_buttons = uigridlayout(bg, [2, 4]);
     rb1 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt1');
     rb2 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt2');
     rb3 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt3');
     rb4 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt4');
     rb5 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt5');
     rb6 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt6');
     rb7 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt7');
     rb8 = uicheckbox(gl2_buttons, ...
         'Text', 'Opt8');
     
     % layout of gl3 (the top panel)
     ax.Layout.Row = [1 6];
     pb1.Layout.Row = 7;
     %layout of gl2 (bottom panel)
     bg.Layout.Row = 5;
     pnl.Layout.Row = [1, 4];
     
     figy.Visible = 'on';
     
     uiwait;
     
     hCheckboxes = findobj(gl2_buttons,'Type','uicheckbox');
     checkboxValues = get(hCheckboxes, 'Value');
     
     close all force
     
     selectedImages = cell2mat(checkboxValues);
     bestIms = ImList(selectedImages); 
     andImage = double(bestIms{1});
     nIms = numel(bestIms);
     for pp = 2:nIms
         andImage = andImage+double(bestIms{pp});
     end
     
     andImage(andImage<nIms)=0;
     andImage(andImage==nIms)=1;
     andImage = logical(andImage);
     gg = regionprops('table', andImage, 'MajorAxisLength');
     nRegions = sum(gg.MajorAxisLength>200);
     finMask = bwareafilt(andImage, nRegions);
     
end

function buttonCallback(~, ~)
    
    uiresume;
end