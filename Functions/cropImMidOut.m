function [cropIm, xRange, yRange] = cropImMidOut(im,varargin)
%function returns a cropped image in the size specified by opt.cropsize+1
%in each dimmension. This function can handle images stored in a 3d matrix

opt.cropsize=[]; % input expects a vector with desired image size in [x y];
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% check if cropsize exceeds im bounds
if opt.cropsize(2) > size(im,1) || opt.cropsize(1) > size(im,2)
    fprintf('cropBounds are larger than the input Image');
    return
end

%%
crop2Range=opt.cropsize/2;

midY=ceil(size(im,1)/2);
midX=ceil(size(im,2)/2);

if ~isempty(opt.cropsize)
    minY=ceil(midY-crop2Range(2));
    maxY=floor(midY+crop2Range(2));
    minX=ceil(midX-crop2Range(1));
    maxX=floor(midX+crop2Range(1));
    im=im(minY:maxY,minX:maxX,:);
end
xRange=[minX maxX];
yRange=[minY maxY];
cropIm=im;
end