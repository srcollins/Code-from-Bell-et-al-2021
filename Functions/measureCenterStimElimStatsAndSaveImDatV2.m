function [imDat, failFlag]= measureCenterStimElimStatsAndSaveImDatV2(...
    correctedRatioIm,grayImg,allOrigMasks,cropsize2,frapMask,allKatIm,allTomIm,allRatioCorr,varargin)
%2020-05-27 GB: This function takes in the fret ratio images pre and post
%application of the fret ratio filter. Additionally, the function takes in
%the corresponding pre and post fret filter masks, as well as the frap mask
%and the cropsize parameters

%The function will crop all input images and save them in imDat data
%structure for future use. Additionally, the function will compute the
%movement statistics that I will use for eliminating cells that move too
%far to be used in the center stim analysis.

%2020-06-04 NOTE!! this function needs editing. I set the function to return a square
%mask that is almost the size of the entire image. When I wrote the
%function we were only concerned with the 10 frames post stim. Now we care
%about the frames before too. Frames pre-stimulation that were replaced
%with the enourmous mask are not caught by this function and these cells
%are not filtered out of the data. I will write a separate function to
%remove these cells from the existing data. 

%2020-06-04 Fixed!
%2020-07-18 V2 removed the images that used the fret ratio as a filter for
%removing pixels.

%2020-07-23 updates required:
%1) update elim stats to handle frap coors as a matrix rather than
%vector-complete!
%2020-07-25 updates: 1)removed andMask tracking, 2) updated frap2edge to
%use original cell mask rather and mask. Also frap2edge now compares each cell
%mask to the corresponding frapmask that has also been shifted. 3) expand
%frap2edge and maskCentroidVect to cover all frames not just 10 post stim.
%euclidDist is for frames after stimulation.

%2020-07-28 Updates: 1) make euclidean distance relative to frame 5 for all
%frames including the before stim frames. Also added andMask tracking back
%per sean's request

%2020-08-05: update: need to also shift the ratioCorrectionFinal images. I
% now have the shifted ratioCorrection images as an input so that they can
% be cropped and saved in imDat.
%% process varargin inputs
opt.numbefore=5;
opt.framerange=[1:length(correctedRatioIm)];


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%% crop Img for saving in imDat and calculate the compute truncation stats;


%crop ratioIm to just barely larger than the cell.
cropFrapIm=cellfun(@(x) cropImMidOut(x,'cropsize',cropsize2),frapMask,'UniformOutput',false);
%frap Images are an array now, need to make the frap coors a matrix.
frapRow=nan(length(cropFrapIm),1); frapCol=nan(length(cropFrapIm),1);
for i=1:length(frapRow)
[frapR, frapC]=find(cropFrapIm{i});
frapRow(i)=frapR; frapCol(i)=frapC;
end
cropfrapInd=[frapCol,frapRow]; %cropFrapInd is now in X,Y format rather than RxC
cropMask=cellfun(@(x) cropImMidOut(x,'cropsize',cropsize2),allOrigMasks,'UniformOutput',false);
cropRatioIm=cellfun(@(x) cropImMidOut(x,'cropsize',cropsize2),correctedRatioIm,'UniformOutput',false);
cropGray=cellfun(@(x) cropImMidOut(x,'cropsize',cropsize2),grayImg,'UniformOutput',false);
cropKat=cellfun(@(x) cropImMidOut(x,'cropsize',cropsize2),allKatIm,'UniformOutput',false);
cropTom=cellfun(@(x) cropImMidOut(x,'cropsize',cropsize2),allTomIm,'UniformOutput',false);
cropRatioCorrection=cellfun(@(x) cropImMidOut(x,'cropsize',cropsize2),allRatioCorr,'UniformOutput',false);
% compute truncation stats; Will use these to filter the data
% later

% calc min dist from frapspot to cell edge

% here I detect if any of the cells have empty masks. I id
% the frame of the emptyMask im and I make it a mask image
% that has the entire image as ones, with a 5 pixel buffer
% as 0s on the margins of the image. If the numbefore im is
% missing, the 3 metrics will have huge values. If one of
% the frame afters is missing, that one frame will have
% huge values. I don't think this will be common, but I did
% catch an error for a movie where the first frame was
% emptyMask.
emptyStructs=cellfun(@(x) isempty(regionprops(x,'Area')),cropMask);
emptyInd=find(emptyStructs);
if ~isempty(emptyInd)
    for id=1:length(emptyInd)
        cropMask{emptyInd(id)}(5:size(cropMask{emptyInd(id)})-5,5:size(cropMask{emptyInd(id)})-5)=1;
    end
end
cropMask=cellfun(@(x) returnMaskWithLargestArea(x),cropMask,'UniformOutput',false);
testAreas=cellfun(@(x) regionprops(x,'Area'),cropMask);
testAreaVect=arrayfun(@(x) x.Area,testAreas);
if all(testAreaVect(opt.framerange)>1000) && all(testAreaVect(opt.framerange)<100000) &&  sum(testAreaVect(opt.framerange)>1000)==length(opt.framerange)    %skips cells with mask areas<1000;
    %compute the distance between the frap spot and the closest point on
    %teh cell edge.
    edgeMask=cell(1,length(cropMask));
    frap2edge=nan(1,length(cropMask));
    for j=1:length(cropMask)
        edgeMask{j}=bwperim(cropMask{j});
        [~,frap2edge(j)]=getClosestPointWithinMask(edgeMask{j},cropfrapInd(j,:));
    end
    frap2edge=round(frap2edge);

    %andMask
    andMask=cellfun(@(x) x&cropMask{opt.numbefore},cropMask,'UniformOutput',false);
    andMask=cellfun(@(x) returnMaskWithLargestArea(x),andMask,'UniformOutput',false);
    andArea=zeros(1,length(andMask));
    for i=1:length(andMask)
        andAreaTemp=regionprops(andMask{i},'Area');
        if isempty(andAreaTemp)
            andArea(i)=0;
        else
            andArea(i)=andAreaTemp.Area;
        end
    end
    andMaskChange=abs(round(((andArea(opt.numbefore)-andArea)/andArea(opt.numbefore))*100));
    
    % calc euclidean distance for the difference between cell mask centroids.
    maskCentroid=cellfun(@(x) regionprops(x,'Centroid'),cropMask);
    maskCentroidVect=nan(length(cropMask),2);
    for f=1:length(maskCentroid)
        maskCentroidVect(f,:)=maskCentroid(f).Centroid;
    end
    maskCentBefore=maskCentroidVect(opt.numbefore,:);
    euclidDist=sqrt((maskCentroidVect(:,1)-maskCentBefore(1,1)).^2 + (maskCentroidVect(:,2)-maskCentBefore(1,2)).^2);
    
    
    
    %save the processed mask and fret images for later
    imDat.cropFrapIm=cropFrapIm;
    imDat.cropfrapInd=cropfrapInd; %XY format
    imDat.cropMask=cropMask;
    imDat.cropFret=cropRatioIm;
    imDat.cropGray=cropGray;
    imDat.cropKat=cropKat;
    imDat.cropTom=cropTom;
    imDat.cropRatioCorrection=cropRatioCorrection;
    imDat.frap2edge=frap2edge;
    imDat.andMask=andMask;
    imDat.andArea=andArea;
    imDat.andMaskChange=andMaskChange;
    imDat.cellCentroids=maskCentroidVect;
    imDat.euclidDist=euclidDist;
    failFlag=false;
else
    failFlag=true;
    imDat=[];
    
end % if sum(testAreaVect(opt.numbefore:opt.numbefore+10)>1000)==11
end