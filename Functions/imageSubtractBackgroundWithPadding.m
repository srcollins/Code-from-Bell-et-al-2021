function [imOut,bg] = imageSubtractBackgroundWithPadding(im,blockSize,padSize,ptile)
%imageSubtractBackgroundWithPadding subtracts the background from image
%block by block, adding an extra region around each block to get more pixel
%data.

if nargin<2
    blockSize=80;
end
if nargin<3
    padSize=60;
end
if nargin<4 | ptile==50
    bg=blockproc(im,[blockSize blockSize],@(x) ones(x.blockSize)*median(x.data(x.data>0)),'BorderSize',[padSize padSize],'TrimBorder',0);
else
    bg=blockproc(im,[blockSize blockSize],@(x) ones(x.blockSize)*prctile(x.data(x.data>0),ptile),'BorderSize',[padSize padSize],'TrimBorder',0);
end
bg=imfilter(bg,fspecial('disk',blockSize/2),'symmetric');
%Subtract the background
imOut=max(im-bg,0);


end
