function [maskFinal,cellCoors]=fretGetCellMasks_63xWithErode(imForSeg0,varargin)
%fretGetCellMasks takes as input an image and a set of bounding rectangles
%(output from fretGetInitialObjectsAndBGmask) and computes a more precise
%mask for all cells detected except those touching the image boundary
opt.objectrectangles={[1,1,size(imForSeg0,2),size(imForSeg0,1)]};
opt.minsize=2000;
opt.imcfp=[];
opt.imfret=[];
opt.sharpparam=[4 50];
opt.erodesize=1;
opt.grad=2;
opt.minperct=5;
opt.maxperct=99.9;

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end



maskFinal=false(size(imForSeg0));
for i=1:length(opt.objectrectangles)
    maskFinal=processSingleObjectRectangle(maskFinal,opt.objectrectangles{i},...
        imForSeg0,opt.minsize,opt.sharpparam,opt.erodesize,opt.grad,opt.minperct,...
        opt.maxperct);
end

objects=regionprops(maskFinal,'BoundingBox','Area','Centroid','PixelIdxList','PixelList');
objects=objects(arrayfun(@(x) x.Area>opt.minsize,objects));

cellCoors(:,1)=vect(arrayfun(@(x) x.Centroid(1),objects));
cellCoors(:,2)=vect(arrayfun(@(x) x.Centroid(2),objects));
cellCoors(:,3)=vect(arrayfun(@(x) x.Area,objects));
for i=1:length(objects)
    cellCoors(i,4:7)=objects(i).BoundingBox;
    if ~isempty(opt.imcfp)&& ~isempty(opt.imfret)
        pixInd=objects(i).PixelIdxList;
        cellCoors(i,10)=sum(opt.imfret(pixInd))/sum(opt.imcfp(pixInd));
        cellCoors(i,11)=median(opt.imfret(pixInd)./opt.imcfp(pixInd));
        
        %CFP and FRET intensities
        cellCoors(i,12)=mean(opt.imcfp(pixInd));
        cellCoors(i,13)=mean(opt.imfret(pixInd));
        
        %Smoothed values for estimating polarization direction
        cfpSmooth=nan(size(opt.imfret)); cfpSmooth(pixInd)=opt.imcfp(pixInd);
        cfpSmooth=ndnanfilter(cfpSmooth,fspecial('gaussian',7,2),7);
        fretSmooth=nan(size(opt.imfret)); fretSmooth(pixInd)=opt.imfret(pixInd);
        fretSmooth=ndnanfilter(fretSmooth,fspecial('gaussian',7,2),7);
        ratioSmooth=fretSmooth./cfpSmooth;
        
        %Weighted centroid:
        cellCoors(i,8)=sum(objects(i).PixelList(:,1).*ratioSmooth(pixInd))/sum(ratioSmooth(pixInd));
        cellCoors(i,9)=sum(objects(i).PixelList(:,2).*ratioSmooth(pixInd))/sum(ratioSmooth(pixInd));
        
        %Measures of polarization
        polarizationVector=[cellCoors(i,8)-cellCoors(i,1) cellCoors(i,9)-cellCoors(i,2)];
        polarizationVector=polarizationVector/norm(polarizationVector);
        relCoord=objects(i).PixelList - repmat(cellCoors(i,1:2),[length(pixInd) 1]);
        projectedCoord=relCoord*polarizationVector(:);
        cellCoors(i,14)=myNanCorrcoef(projectedCoord,ratioSmooth(pixInd));
        ptile75=prctile(projectedCoord,75);
        ptile25=prctile(projectedCoord,25);
        pixRear=pixInd(projectedCoord<=ptile25);
        pixFront=pixInd(projectedCoord>=ptile75);
        cellCoors(i,15)=sum(opt.imfret(pixRear))/sum(opt.imcfp(pixRear));
        cellCoors(i,16)=sum(opt.imfret(pixFront))/sum(opt.imcfp(pixFront));
    end
end

%------------------------------------------------------------------------------
function maskFinal=processSingleObjectRectangle(maskFinal,objRect,imForSeg0,...
    minSize,sharpParam,erodeSize,gRad,percentMin,percentMax)


width=fretGetObjectRectanglePaddingSize;

cx1=floor(max(1,objRect(1)));
cy1=floor(max(1,objRect(2)));
cx2=ceil(min(size(imForSeg0,2),sum(objRect([1 3]))));
cy2=ceil(min(size(imForSeg0,1),sum(objRect([2 4]))));

%Smooth and sharpen the image before segmenting
imForSeg=imForSeg0(cy1:cy2,cx1:cx2);
imForSeg=imfilter(imForSeg,fspecial('gaussian',5,gRad),'symmetric');
imSumNorm=mat2gray(imForSeg, [prctile(imForSeg(:),percentMin) prctile(imForSeg(:),percentMax)]);
imSumNorm=imsharpenSean(imSumNorm,sharpParam(1),sharpParam(2)); %Try enhancing the edges


%Segment the image

threshSeg=graythresh(imSumNorm(imSumNorm>0 & imSumNorm<1));
maskThisObj=imSumNorm>threshSeg;

maskThisObj=bwareaopen(maskThisObj,25,4);
maskThisObj=imdilate(maskThisObj,strel('disk',1));
maskThisObj=imerode(maskThisObj,strel('disk',erodeSize));
%maskThisObj=imclearborder(maskThisObj);
maskThisObj=imfill(maskThisObj,'holes'); %plotCDFs({imSumNorm(:),[threshSeg]});
maskThisObj=bwareaopen(maskThisObj,minSize);

objFound=regionprops(maskThisObj,'BoundingBox');

if length(objFound)>1
    maskFinal(cy1:cy2,cx1:cx2)=maskFinal(cy1:cy2,cx1:cx2) | maskThisObj;
else
    maskFinal(cy1:cy2,cx1:cx2)=maskFinal(cy1:cy2,cx1:cx2) | maskThisObj;
end
