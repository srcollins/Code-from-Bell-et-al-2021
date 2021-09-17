function [maskFinal,correctedRatioIm, grayImg,acceptorIm,donorIm,minAreaThreshPass]= processFretRatioImTIRF(cropIm2,centerMask,ratioCorrectionFinal,varargin)
%2020-05-27 GB: This function runs our normal 60x tirf fret ratio image
%processing. The script takes up lots of space, so I wanted to convert this
%to a function. 

%Inputs: 1) cropIm2: this input is a 3d matirix that contains the donor and
%acceptor fret images that have already been processed for the camera
%darkstate, halfchip and donut corrections. Additionally, these images
%should be registered and cropped.

%2)centerMask: this is a logical mask of the frapspot that has been
%radially dilated to by 50.

%3)Varargin: so far these are just the extra inputs to the masking function. 

%OutPuts: The funtion outputs the following:
% 1)the final Mask
% 2)the corrected ratio image
% 3) gray image that is made by summing the two fluorescent channels


%% process varargin inputs
opt.minsizebg=1200; % for initial bg subtraction
opt.minsize=2000; % for fretGetCellMasks63xErode.
opt.sharpparam=[4 25];
opt.erodesize=1;
opt.grad=1;


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end


%% Process the fret images

% background subtraction and Image mask processing
imSum=sum(cropIm2,3);
[maskBGinitial, objRect ,minAreaThreshPass]=fretGetInitialObjectsAndBGmask4GB(imSum,opt.minsizebg); %checked
if minAreaThreshPass
imForSeg0=imageSubtractBackgroundWithObjectMaskingAndPadding(imSum,150,100,maskBGinitial); %checked
[maskFinal,~]=fretGetCellMasks_63xWithErode...
    (imForSeg0,'minsize',opt.minsize,'imcfp', imSum,'imfret',imSum,...
    'sharpparam',opt.sharpparam,'grad',opt.grad,'erodesize',opt.erodesize,...
    'objectrectangles',objRect); %sharparam is different from default. find out what things in above function are doing
maskFinal=bwareaopen(maskFinal,opt.minsize);
CFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(cropIm2(:,:,2),150,100,maskBGinitial);
YFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(cropIm2(:,:,1),150,100,maskBGinitial);

%restrict cell mask to just the cell that overlaps with the dilated
%frapMask
tempMask=bwlabel(maskFinal);
%tempMask=bwlabel(imdilate(maskFinal, strel('disk',5)),8); % I am not sure why we were dilating this mask
tempMask=tempMask==median(tempMask(tempMask & centerMask));
maskFinal=maskFinal&tempMask;



%smoothing
YFPsmooth=YFPFinal;
CFPsmooth=CFPFinal;
YFPsmooth(~maskFinal)=nan;% mask background from cell
CFPsmooth(~maskFinal)=nan;
YFPsmooth=ndnanfilter(YFPsmooth,fspecial('gaussian',7,2),7);
CFPsmooth=ndnanfilter(CFPsmooth,fspecial('gaussian',7,2),7);
YFPsmooth(~maskFinal)=nan; %mask cell from background
CFPsmooth(~maskFinal)=nan;
% ratio image
ratioImage=YFPsmooth./CFPsmooth;
correctedRatioIm=ratioImage./ratioCorrectionFinal;
grayImg=imSum;
acceptorIm=YFPsmooth;
donorIm=CFPsmooth;


else
    fprintf('processFretRatioImTIRF: MinCellSizeThresh Failure\n');
    maskFinal=[];
    correctedRatioIm=[];
    grayImg=imSum;
    acceptorIm=[];
    donorIm=[];
    return
end