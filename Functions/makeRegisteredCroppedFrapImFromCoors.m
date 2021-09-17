function [frapMask1PX,frapInd,frapMask1PXrt, frapIndrt] = makeRegisteredCroppedFrapImFromCoors(frapPos, p, varargin)

%The original frap spot coordinates are reported in [x y];
% The final frap coors (regFrapInd) are exported in [X Y];
opt.deltaspot=[0 0];
opt.cropsize=[];

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% Build frapMask1PX and frapInd
regFrapInd=zeros(1,2);
frapMask=zeros(1024,1024);
frapMask(frapPos(2),frapPos(1))=1; % remember that mat to image requires (Y X) org
frapMask(:,:,1)=frapMask;
frapMask(:,:,2)=frapMask(:,:,1);
regMask=registerImagesFromQuadFit(frapMask,p);
regFrapMask=regMask(:,:,1);
regFrapMask=cropImMidOut(regFrapMask,'cropsize', opt.cropsize);
[~,frapCenter]=max(vect(regFrapMask(:,:,1)));

[regFrapInd(2),regFrapInd(1)]=ind2sub(size(regFrapMask(:,:,1)),frapCenter);
frapInd=regFrapInd;
frapMask1PX=logical(regFrapMask);
%% build frapMask1PXrt and frapIndrt
frapPos1=frapPos+opt.deltaspot;
regFrapInd=zeros(1,2);
frapMask=zeros(1024,1024);
frapMask(frapPos1(2),frapPos1(1))=1; % remember that mat to image requires (Y X) org
frapMask(:,:,1)=frapMask;
frapMask(:,:,2)=frapMask(:,:,1);
regMask=registerImagesFromQuadFit(frapMask,p);
regFrapMask=regMask(:,:,1);
regFrapMask=cropImMidOut(regFrapMask,'cropsize', opt.cropsize);
[~,frapCenter]=max(vect(regFrapMask(:,:,1)));

[regFrapInd(2),regFrapInd(1)]=ind2sub(size(regFrapMask(:,:,1)),frapCenter);
frapIndrt=regFrapInd;
frapMask1PXrt=logical(regFrapMask);
end