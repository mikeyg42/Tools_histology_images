function [theta_clockwise, flip_LR] = RGB_reorientationGUI(imToRotate,fixedImage)
% syntax = [theta_clockwise, flip_LR] = RGB_reorientationGUI(imToRotate,fixedImage)

%% displays a template (ie the FIXED image) to the left in one axis, adjacent to it 
%% a vertically stacked montage of the MOVING image, rotated;
%% 0º, 90ª, 180º, or 270º (ordered top to bttm), with 4 accompanying buttons to select best.
%% There is also the option to take the mirror image of the moving image 
% NOTE: it will APPEAR as though the FIXED image is being inverted into
% its mirror image... this is an ILLUSION! nothing more... no need to panic!
% Please rest easy knowing any transformations will only be applied to the MOVING image. 

if size(imToRotate, 3)~=3 
    disp('input is supposed to be RGB.... stacking 4 copies of the grayscale input atop each other ');
end

% we will first generate all 4 options 
sz = size(imToRotate, 1:2);
BigSide = max(sz(1), sz(2)); % [size does [rows, column]... [1, 3] is 1x3 is [1, 1, 1]
if sz(1)<sz(2) %we here INSIST nrows>ncols ie height > width
   imToRotate = rot90(imToRotate, -1); % "...make a burger.... not a hotdog."
   fixedImage = rot90(fixedImage, -1);
end
% redefine sz to reflect the change in height
sz = size(imToRotate, 1:2);

% Here we construct a 4-by-1 (stacked vertically) montage of the 4
% potential rotations that the moving image could make to match with the
% template. 
padamt = floor((BigSide-min(sz(1),sz(2)))/2);
imBiggerSide = BigSide;
imReallyLongSide = BigSide*4;

%preallocate with ones..must be ONES function, otherwise you're
%background's high contrastt will become the fixation of all computer vision... the talk of the town! 
if ndims(imToRotate) == 3
    imag0 = ones(imBiggerSide, imBiggerSide, 3, 'double'); % is also rows/columns, i/e rows x cols
    imag0(1:BigSide, padamt:padamt+sz(2)-1, :) = imToRotate(:,:,:);
    rotatedImageTiled = ones(imReallyLongSide, imBiggerSide, 3, 'double');
else
    imag0 = ones(imBiggerSide, imBiggerSide, 1, 'double');
    imag0(1:BigSide, padamt:padamt+sz(2)-1, 1) = imToRotate(:,:, 1);
    rotatedImageTiled = ones(imReallyLongSide, imBiggerSide, 1, 'double');
end

% stack 4 squares into one tall rectangle image


rotatedImageTiled(1:imBiggerSide,:, :) = imag0;
rotatedImageTiled(imBiggerSide+1:imBiggerSide*2, :, :) = permute(imag0(imBiggerSide:-1:1, :, :), [2 1 3]); %90degrees CW
rotatedImageTiled(2*imBiggerSide+1:imBiggerSide*3,:, :) = imag0(imBiggerSide:-1:1,imBiggerSide:-1:1, :); %180degrees
rotatedImageTiled(3*imBiggerSide+1:imBiggerSide*4,:, :) = permute(imag0(:,imBiggerSide:-1:1, :), [2 1 3]); %90degrees CCW

% huge downsample for efficiency of displaying the GUI!!!
rotatedImageTiled = imresize(rotatedImageTiled, 0.1);

%  OPEN the GUI
figRotation = uifigure(1);

nameArray = {'WindowState','HandleVisibility','Name','Visible','NumberTitle', 'IntegerHandle'};
valueArray = {'normal', 'on', 'Reorientation_GUI','off', 'on', 'on'};
set(figRotation, nameArray, valueArray);

gl = uigridlayout(figRotation, [2, 4]);
gl.ColumnWidth = {'5x','5x','6x', '3x'};
gl.RowHeight = {'1x', 50};
axBig = uiaxes('Parent', gl);
axBig.Layout.Column = [1, 2];
axBig.Layout.Row = 1;

axSmall = uiaxes('Parent', gl);
axSmall.Layout.Column = 3; axSmall.Layout.Row =[1, 2];

hold on
imshow(rotatedImageTiled,'Parent',axSmall, 'Border', 'tight');
imshow(fixedImage,'Parent',axBig, 'Border', 'tight');
hold off

pan = uipanel(gl);
pan.Layout.Column = 4;
pan.Layout.Row = [1, 2];
glINSIDEpan= uigridlayout(pan, [4, 1]);
choices = {'Option1', 'Option2', 'Option3', 'Option4'};
rotated0 = uibutton(glINSIDEpan, 'Text', choices{1},'ButtonPushedFcn',@rot0_Callback);
rotated1 = uibutton(glINSIDEpan, 'Text', choices{2},'ButtonPushedFcn',@rot1_Callback);
rotated2 = uibutton(glINSIDEpan, 'Text', choices{3},'ButtonPushedFcn',@rot2_Callback);
rotated3 = uibutton(glINSIDEpan, 'Text', choices{4},'ButtonPushedFcn',@rot3_Callback);
rotated0.Layout.Row = 1; rotated1.Layout.Row = 2; rotated2.Layout.Row = 3; rotated3.Layout.Row = 4;

flipButton = uibutton(gl, 'Text', 'Flip L/R (mirror)', 'ButtonPushedFcn', @flipLR_Callback);
flipButton.Layout.Row = 2;
flipButton.Layout.Column = 2;

set(figRotation, 'Visible', 'on');  % the figure needs to be visible to be made always on top

%initialize the appdata 
setappdata(figRotation, 'flag', 0); %flag is for the flip LR functionality
setappdata(figRotation,'myRotationSelection', 0); %initialize the rotation button

drawnow limitrate % updates figures, but defers any callbacks!

uiwait; %uiresume is hidden away in each of the button callbacks!

%if the uiresume that breaks the uiwait is in the callback for the flipLR
%button, then the appdata would have been changed abd not equal 0 anymore
flip_LR = getappdata(figRotation, 'flag');
if flip_LR ~= 0
    set(figRotation, 'Visible', 'off');
    cla(axBig, 'reset'); 
    
    fixedFlipped = fliplr(fixedImage);
    imshow(fixedFlipped, 'Parent', axBig, 'Border', 'tight');
    
    set(flipButton, 'Visible', 'off'); %we dont have functionality here to permit one to change their mind on that choice. 
    set(figRotation, 'Visible', 'on');
    
    uiwait;
end

choiceNum = getappdata(figRotation, 'myRotationSelection');
choiceStr = num2str(choiceNum);

switch choiceStr
    case '1'
        theta_clockwise = 0;
    case '2'
        theta_clockwise = 90;
    case '3'
        theta_clockwise = 180;
    case '4'
        theta_clockwise = 270;
end

close 2: Reorientation_GUI

end


function rot0_Callback(hObject, ~)
 mygridlayout = hObject.Parent;
 mypanel = mygridlayout.Parent;
 gl2 = mypanel.Parent;
 figRotation = gl2.Parent;
setappdata(figRotation, 'myRotationSelection', 1);
    
uiresume;
end

function rot1_Callback(hObject, ~)
 mygridlayout = hObject.Parent;
 mypanel = mygridlayout.Parent;
 gl2 = mypanel.Parent;
 figRotation = gl2.Parent;
 setappdata(figRotation, 'myRotationSelection', 2);
    
 uiresume;
end

function rot2_Callback(hObject, ~)
 mygridlayout = hObject.Parent;
 mypanel = mygridlayout.Parent;
 gl2 = mypanel.Parent;
 figRotation = gl2.Parent;
 setappdata(figRotation, 'myRotationSelection', 3);
    
 uiresume;
end

function rot3_Callback(hObject, ~)
 mygridlayout = hObject.Parent;
 mypanel = mygridlayout.Parent;
 gl2 = mypanel.Parent;
 figRotation = gl2.Parent;
 setappdata(figRotation, 'myRotationSelection', 4);
    
 uiresume;
end

function flipLR_Callback(hObject, ~)
 mygridlayout = hObject.Parent;
 mypanel = mygridlayout.Parent;
 gl2 = mypanel.Parent;
 figRotation = gl2.Parent;

setappdata(figRotation, 'flag',1)

uiresume;
end
