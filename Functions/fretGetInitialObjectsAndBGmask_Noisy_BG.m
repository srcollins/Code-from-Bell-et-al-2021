function [bgMask,objectRectangles,labeledIm] = fretGetInitialObjectsAndBGmask_Noisy_BG(imSum,minCellSize)
%fretGetInitialObjectsAndBG does an initial pass object detection and
%background subtraction to identify objects for finer processing.

if nargin<2
    minCellSize=2000;
end

%maxValue = 5;
imBS = (imSum-min(imSum(:)))/nanmedian(imSum(:));
% trim so that imBS is square
%imBS = squareIm(imBS);
imBS(imBS>prctile(imBS(:),99.9)) = prctile(imBS(:),99); % reduce highest pixels to preserve the threshold

% generate noise image
imSquare = imBS.*imBS;

%Do an initial segmentation of background and foreground
imSquare=imfilter(imSquare,ones(3)/9,'symmetric');
imSquare = imclearborder(imSquare);

bgMask=imSquare<multithresh(imSquare)*0.8;
bgMask = imerode(bgMask,strel('disk',15)); % strategy - over-expand the crude mask to get around splitting events

labeled = bwlabel(~bgMask);
objects=regionprops(labeled,'Area');


areaThresh = [objects.Area]>=minCellSize;
objects = objects(areaThresh);
areaMask = labeled;

for iArea = 1:length(unique(areaMask))-1
    if areaThresh(iArea) ~= 1
        aMask = areaMask == iArea;
        areaMask(aMask) = 0;
    else
    end % if
end % for iArea

bgMask(~areaMask) = 1;

%bgMask = imclearborder(bgMask == 0);
bgMask = imerode(bgMask,strel('disk',20));
labeledIm = bwlabel(~bgMask);
fgMask = ~bgMask;

objFinal = regionprops(labeledIm,'BoundingBox');



% figure;
% imagesc(imSum); hold on;
for i=1:length(objFinal)
    width=fretGetObjectRectanglePaddingSize*1.5;
    objectRectangles{i}=objFinal(i).BoundingBox+[-1*width -1*width +3*width +3*width]; % taller than wider 
%     rectangle('Position',objectRectangles{i});
end

end

