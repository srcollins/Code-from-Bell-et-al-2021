function resOut=createResFromImgDat(masterRoot, varargin)
%2020-10-27 GB: This function will create the res struct array that
%contains important information for sifting through the centerStim data. I
%am creatign this function as I forgot to save the res variables from a few
%experiments in late september and early october 2020.

%% process varargin inputs:
opt.dirkeyword='Fast';
opt.imdatkeyword='imgData.mat';
opt.numbefore=10;
opt.obj=60;
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%% determine the types of stimulation conditions

rawDataDir=getDirectories([masterRoot 'Analyzed Data'],opt.dirkeyword );
for numRawDir=1:length(rawDataDir)
    parentFolder=rawDataDir{numRawDir};
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    datNames=getFilenames(saveRoot,opt.imdatkeyword); % get the subdirectories for each folder in rawData
    datNames=natsortfiles(datNames);
    cellLabelsAll=regexp(datNames,'(?<=PLB-).*(?=Cell)','match');
    cellTypeLabel=unique(cellfun(@(x) x,cellLabelsAll));
    cellTypeLabels{numRawDir}=cellTypeLabel;
end

count=0;
for i=1:length(cellTypeLabels)
    for j=1:length(cellTypeLabels{i})
      count=count+1;  
cellTypeLabelsCat{count}=cellTypeLabels{i}{j};
    end
end

cellTypeLabelAll=unique(cellTypeLabelsCat);

counter=zeros(1,length(cellTypeLabelAll));


%% preallocate resOut

for i=length(cellTypeLabelAll):-1:1
     resOut(i).name=[];
     resOut(i).dist=[];
     resOut(i).ref5DiffExp=[];
     resOut(i).ref5NumExp=[];
     resOut(i).ref2_5DiffExp=[];
     resOut(i).ref2_5NumExp=[];
     resOut(i).ref2DiffExp=[];
     resOut(i).ref2NumExp=[];
     resOut(i).wholeCellMeanFret=[];
     resOut(i).times=[];
     resOut(i).frapTimeAtZero=[];
     resOut(i).andMaskChange=[];
     resOut(i).frap2edge=[];
     resOut(i).cellCentroids=[];
     resOut(i).euclidDist=[];
     resOut(i).parentFolder=[];
     resOut(i).date=[];
     resOut(i).fileName=[];
     resOut(i).numBefore=[];
end


%% load imDats and compute res
rawDataDir=getDirectories([masterRoot 'Analyzed Data'],opt.dirkeyword );
umPerPix=distanceScale(opt.obj);

for numRawDir=1:length(rawDataDir)
    parentFolder=rawDataDir{numRawDir};
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    datNames=getFilenames(saveRoot,opt.imdatkeyword); % get the subdirectories for each folder in rawData
    datNames=natsortfiles(datNames);
    
    for k=1:length(datNames)
        clear imDat;
        load([saveRoot filesep datNames{k}]);
        
        
        cellTypeLabel=regexp(imDat.fileName,'(?<=PLB-).*(?=Cell)','match');
        %generate logical map for the identified stimGroup.
        cellStimLogicals=strcmpi(cellTypeLabelAll,cellTypeLabel);
        stimInd=find(cellStimLogicals);
        counter(stimInd)=counter(stimInd)+1;
        
        % Compute FRET ratio changes as a function of distance from the frap spot
        distVals=0:0.66:12; % this is in um already
        binWidth=(distVals(3)-distVals(2))/2;
        frameNumRange=[2:length(imDat.cropMask)];
        
        % fill in resOut
        
        %ChunkyBins
        resOut(stimInd).wholeCellMeanFret{counter(stimInd)}=imDat.wholeCellMeanFret;
            %image time
            resOut(stimInd).times{counter(stimInd)}=imDat.times.tkTime';
            frapTimeAtZero=imDat.times.tkInSec-imDat.times.frapInSec;
            resOut(i).frapTimeAtZero{counter(stimInd)}=frapTimeAtZero;
            % truncation metrics
            resOut(stimInd).andMaskChange{counter(stimInd)}=imDat.andMaskChange;
            resOut(stimInd).frap2edge{counter(stimInd)}=imDat.frap2edge;
            resOut(stimInd).cellCentroids{counter(stimInd)}=imDat.cellCentroids;
            resOut(stimInd).euclidDist{counter(stimInd)}=imDat.euclidDist;
            %otherDeets
            resOut(stimInd).parentFolder{counter(stimInd)}=parentFolder;
            resOut(stimInd).date{counter(stimInd)}=imDat.date;
            resOut(stimInd).fileName{counter(stimInd)}=imDat.fileName;
            resOut(stimInd).numBefore(counter(stimInd))=opt.numbefore;
        
        
        
        
        
        for d=1:length(frameNumRange)
            fNumBefore=find(frameNumRange==opt.numbefore);
            fCtrl=find(frameNumRange==2);
            
            %ref5 calc
            ref5DistMask=bwdistgeodesic(imDat.cropFret{frameNumRange(fNumBefore)}>0,imDat.cropFrapIm{frameNumRange(fNumBefore)},'quasi-euclidean')*umPerPix;
            ref5RatioDiff=imDat.cropFret{frameNumRange(d)}-imDat.cropFret{frameNumRange(fNumBefore)};
            [ref5YOut,ref5NOut] = evalBySlidingBins(vect(ref5DistMask),vect(ref5RatioDiff),distVals,binWidth,'nanmean');
            
            %ref5 distMask and ref2 images
            ref2_5RDiff=imDat.cropFret{frameNumRange(d)}-imDat.cropFret{frameNumRange(fCtrl)};
            [ref2_5YOut,ref2_5NOut] = evalBySlidingBins(vect(ref5DistMask),vect(ref2_5RDiff),distVals,binWidth,'nanmean');
            
            %ref2 ctrl
            centMask=false(size(imDat.cropFrapIm{frameNumRange(fCtrl)}));
            obj=regionprops(imDat.cropMask{frameNumRange(fCtrl)},'Centroid');
            centMask(round(obj(1).Centroid(2)),round(obj(1).Centroid(1)))=true;
            ref2DistMask=bwdistgeodesic(imDat.cropMask{frameNumRange(fCtrl)}>0,centMask>0)*umPerPix;
            ref2RDiff=imDat.cropFret{frameNumRange(d)}-imDat.cropFret{frameNumRange(fCtrl)};
            [ref2YOut,ref2NOut] = evalBySlidingBins(vect(ref2DistMask),vect(ref2RDiff),distVals,binWidth,'nanmean');
            
            %store all this data in a sturct array
            resOut(stimInd).name{counter(stimInd)}=cellTypeLabel;
            resOut(stimInd).dist(1,:)=distVals;
            %ref5
            resOut(stimInd).ref5DiffExp{d}(counter(stimInd),:)=ref5YOut;
            resOut(stimInd).ref5NumExp{d}(counter(stimInd),:)=ref5NOut;
            %ref2_5
            resOut(stimInd).ref2_5DiffExp{d}(counter(stimInd),:)=ref2_5YOut;
            resOut(stimInd).ref2_5NumExp{d}(counter(stimInd),:)=ref2_5NOut;
            %ref2
            resOut(stimInd).ref2DiffExp{d}(counter(stimInd),:)=ref2YOut;
            resOut(stimInd).ref2NumExp{d}(counter(stimInd),:)=ref2NOut;
            
        end %d=1:length(frameNumRange)
    end % k=1:length fileNames of imDat.
    
end % rawnumdir
end