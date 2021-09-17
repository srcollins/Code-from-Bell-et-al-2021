function [imOut,bg]=imageSubtractBGScaleEmptyWell(im,pixMask,emptyYFP,emptyCFP,varargin)
%2020-02-07. GB modified the imageSubstractBGscratch function to use this
%for PDL analysis. The main changes will be to how the image variables are
%stored (eg now imsYFP(:,:,lenght(frames)) rather than
%im(:,:,length(channels));). I wil also add a section to crop the images to
%a small central region for scaling purposes, however, the function will
%return the output images in the size of im.

%This function is designed to subract the background in 20X PDL assay
%experimetns. The function takes the image to be BGsubtracted, the bg
%mask(pix) and the empty well images for each imageing channel. The function
%will use the pix BG mask to ID bg pixels in the input image and check the
%mean intensity from these pixels. THe same process will be repeated for
%the empty well image. The empty well image will then be scaled to match
%the brightness in the input image. Finally the intensity scaled empty well
%image will be subracted from the input image. 

 
opt.yfp=false;
opt.cfp=false;
opt.cropsize=[250 250];

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% determine if the input image is YFP or CFP
if opt.yfp
    imYFP=im;
    %crop to small sq around img center for scaling-This doesnt work as
    %well as the whole image, it overestimates the scalefactor
%     cYFP=cropImMidOut(imYFP,'cropsize',opt.cropsize);
%     cPix=cropImMidOut(pixMask,'cropsize',opt.cropsize);
    imYFP(~pixMask)=nan; % set non-scratch pixels to nan
    sMedian=floor(nanmedian(imYFP(:))); % take median of pixels in pix.
    %eYFP=cropImMidOut(emptyYFP,'cropsize',opt.cropsize);
    eYFP=emptyYFP;
    eYFP(~pixMask)=nan; % set non-scratch pixels to nan in the empty image
    eMedian=floor(nanmedian(eYFP(:)));% take median of pixels in the scratch.
    sFactor=sMedian/eMedian;
    scaledEYFP=emptyYFP*sFactor; % scale the emptyYFP image for BG subtraction.
    bg=scaledEYFP;
    imSub=im-scaledEYFP;
    %now set all negative pixels to 0 and add abs(max for pixels <0)
%     negPix=imSub<0;
%     minBG=abs(min(imSub(negPix)));
%     %imSub(negPix)=1;
%     imOut=imSub+(minBG+1);
    imOut=imSub;
    
        
elseif opt.cfp
    imCFP=im;
    imCFP(~pixMask)=nan; % set non-scratch pixels to nan
    sMedian=floor(nanmedian(imCFP(:))); % take median of pixels in pix.
    %eYFP=cropImMidOut(emptyYFP,'cropsize',opt.cropsize);
    eCFP=emptyCFP;
    eCFP(~pixMask)=nan; % set non-scratch pixels to nan in the empty image
    eMedian=floor(nanmedian(eCFP(:)));
    sFactor=sMedian/eMedian;
    scaledECFP=emptyCFP*sFactor; % scale the emptyYFP image for BG subtraction.
    bg=scaledECFP;
    imSub=im-scaledECFP;
    %now set all negative pixels to 0 and add abs(max for pixels <0)
%     negPix=imSub<0;
%     minBG=abs(min(imSub(negPix)));
%     %imSub(negPix)=0;
%     imOut=imSub+(minBG+1);
    imOut=imSub;
    
elseif opt.yfp==false && opt.cfp==false
    fprintf('select whether the image is a YFP or CFP image');
    return
end





end