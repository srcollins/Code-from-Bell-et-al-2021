function [xRange, yRange, bBoxes]=getXYRangeForImMontageCrop(ImgIn,imSelector,varargin)
% 2020-04-08GB
%This function is designed to take images from a sequence, track a cell,
% and record the bounding box for the cell. The bounding boxes are then
% used to determine the xy range required for croping images to fit the
% image montage.

%inputs: The function will take an image cell array. I was anticipating
%giving it all of the ratio images or all of the grayscale images. The processed
% fret or gray images will have the bg set to nans. The function assums that nans are bg pix and
% that non nan pix are cell pix. 
%frames for the montage are specified in imSelector. imSelector should be a
%vector that contains the frame numbers to be used in the montage. 


opt.padtopleft=50;
opt.padbotright=25;
opt.pos1=[];

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%%
for i=1:length(imSelector)
  clear cellMask;  
cellMask=ImgIn{imSelector(i)};
cellMask(~isnan(cellMask))=1;
cellMask(isnan(cellMask))=0;

if ~isempty(opt.pos1)
   cellMask=shiftMatrix(cellMask,opt.pos1(imSelector(i),1),opt.pos1(imSelector(i),2));
  cellMask1{i}=shiftMatrix(cellMask,opt.pos1(imSelector(i),1),opt.pos1(imSelector(i),2));
end

cellBox(i)=regionprops(cellMask,'BoundingBox');
cellBox(i).BoundingBox=round(cellBox(i).BoundingBox,-1);

bBoxes{i}=cellBox.BoundingBox;  
end
    xRange(1)=min(arrayfun(@(x) x.BoundingBox(1),cellBox))-opt.padtopleft;
    xRange(2)=(max(arrayfun(@(x) x.BoundingBox(1)+x.BoundingBox(3),cellBox)))+opt.padbotright;
    yRange(1)=min(arrayfun(@(x) x.BoundingBox(2),cellBox))-opt.padtopleft;
    yRange(2)=max(arrayfun(@(x) x.BoundingBox(2)+x.BoundingBox(4),cellBox))+opt.padbotright;
      

end

