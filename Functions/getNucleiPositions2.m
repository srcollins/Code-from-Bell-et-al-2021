function [coors,pix,mask]=getNucleiPositions2(im,varargin)
%takes a grayscale image and finds the centers of the nuclei of my HL60
%cells

%20201218 GB Updates:
% I am using this function for segmenting PLB Cdc42-TomKat cells in the 20x
% PDL assay. The goal is to do single cell tracing. In the original version
% of this function, the image threshold was given as function input. Insead
% I will have the function do the thresholding and masking. I used the guts
% of the processSingleObjectRectangle function from
% fretGetCellMaskst63xErode. 

% 2021-03-03 Updates: currently, im should be a gray scale image. I am
% adding an if loop that will check whether im is logical or not.

%% process varargin inputs
opt.minarea=2;
opt.maxarea=Inf;
opt.minsize=50; % this is for bwareaopen, which removes the objects smaller than this value.
opt.maxsize=450; %this is used for bwareaopenBigObj which removes the largest objects in the mask
opt.sharpparam=[4 15];
opt.erodesize=1;
opt.grad=2;
opt.minperct=5;
opt.maxperct=99.9;
opt.imclearborder=true;

opt.largefilt=0; % max pixel number in object, P input for bwareaopen.
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
% %% Defaults


%% Make the mask

%check that im is a grayscale image that needs masking
if ~islogical(im)
    
    %Smooth and sharpen the image before segmenting
    imForSeg=im;
    imForSeg=imfilter(imForSeg,fspecial('gaussian',5,opt.grad),'symmetric');
    imSumNorm=mat2gray(imForSeg, [prctile(imForSeg(:),opt.minperct) prctile(imForSeg(:),opt.maxperct)]);
    imSumNorm=imsharpenSean(imSumNorm,opt.sharpparam(1),opt.sharpparam(2)); %Try enhancing the edges
    
    %Segment the image
    %threshSeg=graythresh(imSumNorm(imSumNorm>0 & imSumNorm<1));
    thresh=multithresh(imSumNorm(imSumNorm>0 & imSumNorm<1),2);
    thresh=thresh(1);
    maskThisObj=imSumNorm>thresh;
    maskThisObj=bwareaopen(maskThisObj,25,4);
    maskThisObj=imdilate(maskThisObj,strel('disk',1));
    maskThisObj=imerode(maskThisObj,strel('disk',opt.erodesize));
    %maskThisObj=imclearborder(maskThisObj);
    maskThisObj=imfill(maskThisObj,'holes'); %plotCDFs({imSumNorm(:),[threshSeg]});
    maskThisObj=bwareaopen(maskThisObj,opt.minsize);
    mask=bwareaopenBigObj(maskThisObj,opt.maxsize);
    if opt.imclearborder
    mask=imclearborder(mask);
    end
else
    mask=im;
end

%Note: I removed the watershed feature as its hard to get it to work well.
%% Find objects
% note - for a sample image which seemed to work well, the median object
% area was 12 pixels. I could adjust the   to shoot for this median
% size
dl=regionprops(mask,'Centroid','PixelIdxList','Area');%get centroid & pixelidxlist
im=double(im);
num=size(dl,1);
coors=nan(num,6);
if nargout>1
    pix=cell(num,1);
end
for l=1:num
    loc=dl(l).Centroid;
    coors(l,1:2)=loc;
    coors(l,3)=median(log10(double(im(dl(l).PixelIdxList))));
    coors(l,4)=dl(l).Area;
    coors(l,5)=log10(sum(vect(double(im(dl(l).PixelIdxList)))));
    coors(l,6)=log10(max(double(im(dl(l).PixelIdxList))));
    if nargout>1
        pix{l}=dl(l).PixelIdxList;
    end
end

ind=coors(:,4)<=opt.maxarea & coors(:,4)>=opt.minarea;
coors=coors(ind,:);
if nargout>1
    pix=pix(ind);
end
%coors(coors(:,4)>opt.maxarea,:)=[];
%coors(coors(:,4)<opt.minarea,:)=[];     %remove objects that are too small
%plotCDFs({coors(:,4)});