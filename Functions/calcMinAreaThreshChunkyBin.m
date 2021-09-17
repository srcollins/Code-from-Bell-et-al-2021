function areaThresh=calcMinAreaThreshChunkyBin(distMasks,numBefore,binRange)
% 2020-06-17 GB: I am using the chunky bin strategy to look at cell mean
% fret ratio based on the distance from the stimulation site. I would like
% to keep track of the area of each bin and exclude cells with areas that
% are smaller than bin1. This function will calculate bin1 area for you
% using the frame before stimulation as the area source. 


%% function meat
tempBin=distMasks{numBefore}>=binRange{1}(1) & distMasks{numBefore}<binRange{1}(2);
a=regionprops(tempBin,'Area');
if isempty(a)
    areaThresh=nan;
else
    areaThresh=a.Area;
end
end