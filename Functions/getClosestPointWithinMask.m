function [pt,minPixelDist]=getClosestPointWithinMask(mask,pt0)
% 2020-05-09 GB edit: added the minPixelDist output so that the function would
% return the distance for the pixel with the closest euclidean distance 
%to the target point.

%2020-07-25: pt0 is in xy format. The function converts to RxC format

ptMask=zeros(size(mask));
ptMask(pt0(2),pt0(1))=1;
dmap=bwdist(ptMask);
dmap(~mask)=Inf;
minVal=min(dmap(:));
minPixelDist=minVal;
[pt(2),pt(1)]=find(dmap==minVal,1);
