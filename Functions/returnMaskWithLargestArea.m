function bigMask=returnMaskWithLargestArea(maskIm)
% 2020-13-05 GB: This function will take a mask inupt image, find the
% largest masked object, and remove the other smaller objects. The function
% will return a logical mask image that only contains the largest object.

%%
tempMask=maskIm;
objs=regionprops(tempMask, 'Area');
if isempty(objs)
    bigMask=zeros(size(maskIm));
else
    [~, maxInd]=max(arrayfun(@(x) x.Area, objs));
    tempMask=bwlabel(tempMask);
    tempMask(tempMask~=maxInd)=0;
    tempMask(tempMask>0)=1;
    bigMask=tempMask;
end
end