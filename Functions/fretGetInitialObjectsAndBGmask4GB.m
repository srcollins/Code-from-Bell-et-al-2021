function [bgMask,objectRectangles,minAreaThreshPass] = fretGetInitialObjectsAndBGmask4GB(imSum,minCellSize)
%fretGetInitialObjectsAndBG does an initial pass object detection and
%background subtraction to identify objects for finer processing.

%2021-02-26 GB update. I added the minAreaTheshPass output. This will
%report whether a cell mask is detected by region props. I am also running
%into an error when regionprops does not detect a cell mask, leaving the
%regionprops output structure empty. I will make it so that if no mask is
%detected, the minAreaThreshPass variable is set to zero and the
%objectrectangels is set to the size of the input image.

if nargin<2
    minCellSize=750;
end

%Do an initial segmentation of background and foreground
%imSum=imfilter(imSum,ones(3)/9,'symmetric');
imNorm=mat2gray(imSum, [prctile(imSum(:),5) prctile(imSum(:),99.5)]);
threshInitialSeg=graythresh(imNorm(imNorm>0 & imNorm<1));
bgMask0=imNorm<threshInitialSeg;
bgMask=imerode(imdilate(bgMask0,strel('disk',5)),strel('disk',20));
fgMask=imerode(imdilate(~bgMask0,strel('disk',5)),strel('disk',10));
fgMask=bwareaopen(fgMask,minCellSize);

%Identify objects and bounding rectangles
objects=regionprops(fgMask,'BoundingBox','Area');
minAreaThreshPass=~isempty(arrayfun(@(x) x.Area, objects));

if isempty(objects) % no cell is detected, set objrect to size of input im and fail the minAreaThreshPass.
    objects(1).BoundingBox={[1,1,size(imSum,2),size(imSum,1)]};
    objectRectangles{1}=objects(1).BoundingBox;
    minAreaThreshPass=0;
    
else
    % figure;
    % imagesc(imSum); hold on;
    for i=1:length(objects)
        width=fretGetObjectRectanglePaddingSize;
        objectRectangles{i}=objects(i).BoundingBox+[-1*width -1*width 2*width 2*width];
        %     rectangle('Position',objectRectangles{i});
    end
end


end