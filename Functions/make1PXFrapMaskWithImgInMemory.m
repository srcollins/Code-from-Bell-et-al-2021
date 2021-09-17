function [frapMask1PX, frapInd] = make1PXFrapMaskWithImgInMemory(img,  varargin)
% Function inputs: frapImage from memory that has been p-registered. Optional inputs include
% crop Range of final Im (ie [600 600]),and translate spot ([x,y];
% Use name value pairs with names coming from the opt structure.
% function returns a logical mask of a 1 pixel frapspot based on the max
% intensity measurement of the actual frap image. The translate spot parameter will shift the spot by the specified pixel distances
% eg [12, 0] = 12 px to the right
% NOTE this function is designed for fret pair experiments and requires a
%a frap image that has already been aligned with the p-alignment parameter.
% Use makeCompositeFrapIm function to build the aligned frap image.
%% varargin inputs
tempRoot=root;

opt.deltaspot=[0 0];
opt.cropsize=[];

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

if nargin < 2
    fprintf('data root and alignment parameters are required.');
    return
end

    %%
    frapIm=img;
   
    if size(frapIm,3)>1
        frapIm=sum(frapIm,3);
    end
    
    frapTemp=frapIm;
    frapTemp=cropImMidOut(frapTemp,'cropsize', opt.cropsize);
    
    frapCenter=find(frapTemp==max(vect(frapTemp)));
    if frapCenter>1
        frapCenter=frapCenter(1);
    end
    [frapInd(1),frapInd(2)]=ind2sub(size(frapTemp),frapCenter);
    % incorporate the deltaSpot into frapInd
    frapInd(1)=frapInd(1)+opt.deltaspot(2);
    frapInd(2)=frapInd(2)+opt.deltaspot(1);
    tempMask=false(size(frapTemp(:,:,1)));
    tempMask(frapInd(1),frapInd(2))=1;
    frapMask1PX=tempMask;
   
end
    
    