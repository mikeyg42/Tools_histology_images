
function visualizeROIs(ROIcentroidCoords, BWimageOfROIS, imgCurrent, replacement, handy)
%% first we overlay 10 roi's on to image
close all
NUMROI = handy.NUMROI;

fHand =figure;
set(fHand,'doublebuffer','off');

imshow(imgCurrent);
set(gca, 'xlimmode','manual',...
    'ylimmode','manual',...
    'zlimmode','manual',...
    'climmode','manual',...
    'alimmode','manual',...
    'Units','pixels');
refreshdata(fHand)

fHand.Visible = 'off';

xlimitIMG = get(gca,'XLim');
ylimitIMG = get(gca,'YLim');
width = floor(xlimitIMG(2));
height = floor(ylimitIMG(2));

magentaRectangle = cat(3, ones([height width], 'double'), zeros([height width], 'double'), ones([height width], 'double'));
gr = cat(3, zeros([height width], 'double'), ones([height width], 'double'), zeros([height width], 'double'));

%prepare the masks to be applied as the alpha data on the greeen layer for visialzation
%and 
les_masques = imresize(BWimageOfROIS, [height width]);
les_masques = double(les_masques);
les_masques = les_masques.*0.4; 

hold on
grIm = image(gr);
hold off

set(grIm, 'AlphaData', les_masques);
refreshdata(gca);

displacement = round(handy.ROISIZE(1)/4);

for jj = 1:NUMROI % Get centroid coords one at a time
    if  ~isempty(ROIcentroidCoords(jj).points)
        rectumPoints = ROIcentroidCoords(jj).points; %a 4x2 array of verticies.
        
        text(rectumPoints(1,1)+displacement, rectumPoints(1,2)+displacement, num2str(jj), 'FontSize', 14, 'FontWeight', 'Bold');
        
    elseif ~isempty(replacement)
        new_masque = imresize(replacement,  [height width]);
        new_masque = double(new_masque).*0.4;
        
        hold on
        magentaImage = image(magentaRectangle);
        hold off
        set(magentaImage, 'AlphaData', new_masque)
    end
    pause(0.1);
end
end


