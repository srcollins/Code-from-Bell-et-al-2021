function distMask=plbMasks2DistMask(masks,frame)

mask0=masks{frame-1};
mask1=masks{frame};
mask2=masks{frame+1};

edge=mask1-imerode(mask1,strel('disk',1));
%diff1=imopen(mask1-mask0,strel('disk',1));
%diff2=imopen(mask2-mask1,strel('disk',1));
diff1=mask1-mask0;
diff2=mask2-mask1;

pEdge=(edge>0) & (diff1>0);
pEdge2=(edge>0) & imdilate(diff2,strel('disk',1));
frontInitial=pEdge & pEdge2;

objects=regionprops(mask1,'PixelIdxList');
front=false(size(frontInitial));
for i=1:length(objects)
    thisMask=false(size(frontInitial));
    thisMask(objects(i).PixelIdxList)=true;
    thisFront=frontInitial & thisMask;
    obj2=regionprops(thisFront,'Area','PixelIdxList');
    if ~isempty(obj2)
        [~,k]=max([obj2.Area]);
        front(obj2(k(1)).PixelIdxList)=true;
    end
end

distMask=bwdistgeodesic(mask1,front,'quasi-euclidean');
%showImagesWithLinkedAxes({mask1,distMask});
