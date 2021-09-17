%% Script for computing dynamics of FRET sensor data for PDL assay
%Updated 2020-08-06: This script is for processing the multi-stim global
%data that has metadata files. The goal is to generate one data structure
%that has all of the multistim data in it with the newest processing
%method. This includes using the dead cell filter and the empty wells for
%assesing bg fluorescence. Additionally the script will extract the time
%stamps from the meta data and save this info in the data structure.
%Finally, I will create a second data structure that contains the mean
%values for both the fret and the exp timing for modeling. 

%2020-10-17 update: I updated how the BG threshold is calculated for pix.
%Originally the script was hard coding bgpixels as any pixel<250 intensity
%units. Some of my LatA data had higher bg levels due to highish cell
%density. As a result, some of the cells had noisy response curves due to
%an empty backgroud. Instead, I am settign the bg threshold as the bottom
%1.5% of pixels. I will check to make sure that this fixes the noise issue.
%v9 of the script also has movie writing capability although I turned it
%off for processing speed. The 1.5% threshold is working nicely.

%2021-03-19 Used this script with 2020-10-21-Cdc42-AllMultiPulse-v1.mat and
%pDatPDL20201022v1.mat data to re plot the medium power ret, no ret control
%plot from fig 1d of the opsin paper. In the original plot, the time stamps
%from when the camera took the picture were not actually used. This caused
%the duration of frames for the stimulated images to be inaccurate.

%2021-03-20 update GB: Currently, the getImgTimesForCellMetaDataPDLassay
%function is putting the image time stamp in the numBefore ref where the
%last frame before stim will have a time stamp of 0. This means that the
%stim time will be a positive number. Instead the time stamps should be in
%the reference of the first stim time so that the stim(1) is occurring at
%t=0. As a temporary work around, I am computing the stimRef times in
%pDatPDL. I added new varargin inputs to buildTimeScaleBarsForPDLPlots to
%allow the user to input the fieldnames of the in ref image times and in
%ref stim times. Sadly I had to remake the PDL figures in fig 4 using these
%updates.

%2021-03-24 Plotted two pulse data for sfig2. 
%2021-03-25 Loaded all Multipulse data. Plotted PAKinhibitor at medium
%power.
%% set roots
masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\GlobalStim-Long\';
%masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\';
%masterRoot = 'C:\Users\George\Documents\temp Data Folder\'; 
%% camera and donut corrections
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat',...
    'donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');

%ratioCorrectionImage
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\2D ratio correction-TK20xMaster.mat');

% load empty BG images
load(['D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\'...
    'global multiStim intensity titration\Analyzed Data\2020-02-04_PLB_Cdc42KO-PakKo-multistim_MedPwr'...
    filesep 'emptyP384WellsTK.mat'],'emptyYFP', 'emptyCFP');
%% define the size of rawDataDir
rawDataDir=getDirectories([masterRoot 'Raw Data'], 'Long');
%rawDataDir=rawDataDir(end);

%% define cLabel, a structure array that contains label info for each cell
% the struct will contain the cell treatment, number of stims and light
% stim power data. The structure is defined outside of the primary data
% processing loop so that all image folders can be identified across
% multiple data directories. The total number of cells will also be used to
% preallocate the dat variable, which is a struct array that contains the
% processed image data.


cLabel=[];
for numRawDir=1:length(rawDataDir)
    clear dataSubdir;
    parentFolder=rawDataDir{numRawDir};
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    dataSubdir=getSubdirectories(rawDataroot,'SC147'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the data
    
    cellLine='PLB';
     constructLabels='SC344-SC147';
    cTemp=getUniqueLabelsPDLAssay(dataSubdir,cellLine,'constructLabels',constructLabels,'keepplasmidlabels',true);
    cLabel=[cLabel cTemp];
    
end
%% adjust cLabels
% sStr={'Cdc42KO','^C10$','^Ctrl$','^LatA$','^LatA-Ret$','Inh','LatA-C10','^LatA-Ctrl$'};
% repStr={'C10-Ret','C10-Ret','Ctrl-Ret','LatA-Ctrl-Ret','LatA-Ctrl-Ret','PAK1Inhib-Ret','LatA-C10-Ret','LatA-Ctrl-Ret'};

%adj for long stim exp
sStr={'SC344-SC147$','SC344-SC147-LatA'};
repStr={'Ctrl-Ret','LatA-Ctrl-Ret'};
cLabelOut=cLabelAdjPDLassay(cLabel,sStr,repStr);

    tCond1= unique(arrayfun(@(x) x.tCond,cLabelOut,'UniformOutput',false));
stimGout= unique(arrayfun(@(x) x.stimGroup,cLabelOut,'UniformOutput',false));


%% preallocate dat
for i=length(cLabel):-1:1
    dat(i).pix_meanA=[];
    dat(i).fileName=[];
    dat(i).meanC=[];
    dat(i).meanY=[];
    dat(i).date=[];
    dat(i).time=[];
    dat(i).tCond=[];
    dat(i).numStim=[];
    dat(i).pwr=[];
    dat(i).dataSubdir=[];
    dat(i).stimGroup=[]; % =[tCond-pwr-numStim];
end
%% compute alignment and ratio correction data for each experimental day


countSubdir=0; % counter for total number of subdirectories read, allows me to build a big struct for 3 days of exp
TFtomKat=true;

imgName1='TomKat-1';
imgName2='TomKat-2';
cropsize=[951 921];
for numRawDir=1:length(rawDataDir)
    clear dataSubdir;
    parentFolder=rawDataDir{numRawDir};
    mkdir([masterRoot 'Scripts' filesep parentFolder]); % makes a folder for the exp day in the scripts folder
    mkdir([masterRoot 'Analyzed Data' filesep parentFolder]);% makes a folder for the exp day in the Analyzed Data folder
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    dataSubdir=getSubdirectories(rawDataroot,'SC344'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the data
    
    
    if exist([scriptRoot filesep 'alignment parameters pX pY.mat'])~=2;
        if TFtomKat
            p0=buildP0TomKat;
        else
            p0=buildP0TomKat('x1',0,'y1',0);
        end
        [p, imgF] = collectImgAndAlignWithPQuadFit(rawDataroot, dataSubdir(1:5), {imgName1,imgName2}, 'p0', p0, 'binsize',1);
        showImagesMergeChannels(imgF(:,:,1),imgF(:,:,2));
        save([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
    else
        % if P is defined, load p parameters.
        load([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
    end % if e if exist([scriptRoot filesep 'alignment parameters pX pY.mat'])==0;
    
    
    
    % align and crop empty well images for bg subtraction
    [regEmptyYFP, regEmptyCFP]=registerEmptyWellImgScratchAssay(emptyYFP,emptyCFP,p);
    cropEyfp=cropImMidOut(regEmptyYFP,'cropsize', cropsize);
    cropEcfp=cropImMidOut(regEmptyCFP,'cropsize', cropsize);
    
    %build metadata file
    metaFrame = parseFrameMetaData(rawDataroot);
    % determine the numbefore img
    numBefore = findNumImgBeforeDapiStimPDLv2(rawDataroot,dataSubdir,{imgName1,imgName2});
    % Fixed bg pixel method FRET ratio analysis
    for subNum=1:length(dataSubdir)
        
        %      for subNum=keepInd{numRawDir}
        clear imTimePre imTimePost stimTimePre stimTimePost;
        countSubdir=countSubdir+1;
        fprintf('i=%i: ',subNum);
        folder=[rawDataroot filesep dataSubdir{subNum}];
        filenames1=getFilenames(folder,imgName1);
        filenames2=getFilenames(folder,imgName2);
        
        
        
        % get time stamps for the current set of fileNames
        if isempty(metaFrame)
            ts=buildTSforPDLdataPreMetaFiles(folder, filenames1);
        else
            ts = getImgTimesForCellMetaDataPDLassay(metaFrame, filenames1,numBefore);
        end
        
        cropsize=[951 921];
        imsY=cell(1,length(filenames1));
        imsC=cell(1,length(filenames2));
%         grayIm=cell(1,length(filenames1));
%         ratioImMovie=cell(1,length(filenames1));
        tic;
        parfor frameNum=1:length(filenames1)
            %fprintf('%i ',frameNum);
            %clear im
            
            im1=double(imread([folder filesep filenames1{frameNum}]));
            im1=(im1- medianDark1).*halfCorrection1;
            im2=double(imread([folder filesep filenames2{frameNum}]));
            im2=(im2-medianDark2).*halfCorrection2;
            im = singleIms2stack(im1,im2);
            regIm=registerImagesFromQuadFit(im,p);
            cropIm=cropImMidOut(regIm,'cropsize', cropsize);
            imsY{frameNum}=cropIm(:,:,1);
            imsC{frameNum}=cropIm(:,:,2);
           % grayIm{frameNum}=imsY{frameNum}+imsC{frameNum};
            
        end
        
        imY=[]; imC=[];
        imY=cat(3,imsY{:});
        imC=cat(3,imsC{:});
        
        
        imAllY=max(imY,[],3);
        imAllC=max(imC,[],3);
        imAllY0=min(imY,[],3);
        imAllC0=min(imC,[],3);
        tempPix=abs(imAllY+imAllC-imAllY0-imAllC0);
        bgThresh=prctile(tempPix,1.5,'all');
        %pix=abs(imAllY+imAllC-imAllY0-imAllC0)<250; %bg == pixels where the max and min are <600 intensity units apart.
        pix=tempPix<bgThresh;
        
        % define the dead cells here use only the frames before stim to
        % find the dead cells and then set to nan
        %sum first 10 frames for doing the ratiocorr on the dead cell filter.
        yDead=nan(size(imAllY,1),size(imAllY,2),10);
        cDead=nan(size(imAllY,1),size(imAllY,2),10);
        for d=1:10
            tempCFP=imC(:,:,d); tempYFP=imY(:,:,d);
            [bgCFP,~]=imageSubtractBGScaleEmptyWell(tempCFP,pix,cropEyfp,cropEcfp,'cfp',true);
            [bgYFP,~]=imageSubtractBGScaleEmptyWell(tempYFP,pix,cropEyfp,cropEcfp,'yfp',true);
            yDead(:,:,d)=bgYFP;
            cDead(:,:,d)=bgCFP;
        end
        
        sumYDead=sum(yDead,3);
        sumCDead=sum(cDead,3);
        tempRatio=sumYDead./sumCDead;
        deadRatio=tempRatio./ratioCorrectionFinal;
        %density_scatter_heatmap(sumYDead+sumCDead,deadRatio,0:100:7e4,0:0.01:3);
        deadThresh=0.8;     %empirically determine this value by looking at the density scatter heat map
        % deadMask=deadRatio<deadThresh & sumCDead >1500;
        %deadMask=deadRatio<deadThresh;
        
        
        
        
        for j=1:length(filenames2)
            tempCFP=imC(:,:,j); tempYFP=imY(:,:,j);
            [bgCFP,~]=imageSubtractBGScaleEmptyWell(tempCFP,pix,cropEyfp,cropEcfp,'cfp',true);
            [bgYFP,~]=imageSubtractBGScaleEmptyWell(tempYFP,pix,cropEyfp,cropEcfp,'yfp',true);
            
            
            %dont take the mean of the ratio image as small and large ratios
            %will bias the result. Instead take the mean of both channels
            %and then ratio.
            corrYFP=bgYFP./ratioCorrectionFinal;
            corrCFP=bgCFP;
            %define dead pixels
            tempRatio=corrYFP./corrCFP;
            
            %            %check frequency distributions for notebook
            %            deadMask=deadRatio<deadThresh;
            %            bigDeadMask=imdilate(deadMask,strel('disk',5));
            %            plotFrequencyDistributions({tempRatio,tempRatio(deadMask),tempRatio(~deadMask),tempRatio(bigDeadMask), tempRatio(~bigDeadMask)},0:0.01:3)
            %             title('C10 Cell')
            %            ylabel('Frequency Density')
            %            xlabel('Fret Ratio')
            
            badMask=tempRatio<deadThresh | tempYFP>6e4 | tempCFP>6e4 | corrYFP<100 | corrCFP<100;
            % plotFrequencyDistributions({tempRatio,tempRatio(badMask),tempRatio(~badMask)},0:0.01:3)
            corrYFP(badMask)=nan;
            corrCFP(badMask)=nan;
%             ratioImMovie{j}=corrYFP./corrCFP;
%              speckleMask=removeSpeckles(corrYFP+corrCFP,'areathresh',50,'pixperbin',20,'intensitythresh',1000);
%             speckleMask=logical(speckleMask);
%              ratioImMovie{j}(speckleMask)=nan;
            
            dat(countSubdir).pix_meanA(j)=nanmean(vect(corrYFP))/nanmean(vect(corrCFP)); %(mean(vect(temp1))-mean(temp1(pix)))/(mean(vect(temp2))-mean(temp2(pix))); this was the older version
            dat(countSubdir).fileName=dataSubdir{subNum};
            dat(countSubdir).meanC(j)=nanmean(vect(corrCFP));
            dat(countSubdir).meanY(j)=nanmean(vect(corrYFP));
            dat(countSubdir).date=regexp(rawDataDir{numRawDir}, '^[0-9\-]+', 'match');
            dat(countSubdir).time=ts;
            dat(countSubdir).tCond=cLabelOut(countSubdir).tCond;
            dat(countSubdir).numStim=cLabelOut(countSubdir).numStim;
            dat(countSubdir).pwr=cLabelOut(countSubdir).pwr;
            dat(countSubdir).dataSubdir=cLabelOut(countSubdir).dataSubdir;
            dat(countSubdir).stimGroup= sprintf('%s-%s-%s',cLabelOut(countSubdir).tCond,cLabelOut(countSubdir).pwr,cLabelOut(countSubdir).numStim);% =[tCond-pwr-numStim];
        end
        
        
%         moviePath=[saveRoot filesep dataSubdir{subNum} '.avi'];
%         writeFRETMoviesFromImageSequences(ratioImMovie,grayIm,moviePath,...
%             'imframerate',ts.frameRate,'objective',20,'boundsfret',[.8 1.2],'scalebarlength',75);
        toc;
        fprintf('\n');
        
    end % subNum=1:length(dataSubdir) ratioVals
    
end % numRawDir
%% 2020-10-21 GB: fix stim time error for files that had timeData as the meta
%data (2019-09-13 and 2020-09-26)

dateInd=[];

for i=1:length(dat)
    dateInd(i)=boolRegExp(dat(i).date,'2019-09');
end

fixInd=1:length(dat);
fixInd=fixInd(logical(dateInd));

for j=1:length(fixInd)
    if ~isempty(dat(fixInd(j)).time.frapTime)
       tNumBefore=abs(dat(fixInd(j)).time.tkTime(1));
       dat(fixInd(j)).time.frapTime= dat(fixInd(j)).time.frapTime-tNumBefore;
    end
end
%%
%tempSaveRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\Analyzed Data\2020-08-07-all-multistimDats';
%tempSaveRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\Analyzed Data\2020-10-12_multiStimDat';
fileName='2021-03-24-Cdc42-AllLong-v1.mat';

save([saveRoot filesep fileName],'dat','-v7.3');

%% save a copy of this script in the scriptRoot

Version='1';
FileNameAndLocation=matlab.desktop.editor.getActiveFilename;
%get just the script name. want to save in scriptRoot location
bslashInds=regexp(FileNameAndLocation,'\\');
fName=FileNameAndLocation(bslashInds(end)+1:end);
fName(end-1:end)=[];
newbackup=sprintf('%s\\%s-backup%s.m',scriptRoot,fName,Version);
A = exist(newbackup,'file')
if (A~=0)
warning('Backup already exists for the current version')
else
    copyfile(FileNameAndLocation,newbackup);
end
%% store dat as a new variable so that I can load the old dat to 
%concatenate the new data into the total dat structure.
oldestDat=dat;
clear dat

%% adjust cLabels
sStr={'Cdc42KO','^C10$','^Ctrl$','^LatA$','^LatA-Ret$','Inh','LatA-C10','^LatA-Ctrl$'};
repStr={'C10-Ret','C10-Ret','Ctrl-Ret','LatA-Ctrl-Ret','LatA-Ctrl-Ret','PAK1Inhib-Ret','LatA-C10-Ret','LatA-Ctrl-Ret'};
cLabelOut=cLabelAdjPDLassay(cLabel,sStr,repStr);

    tCond1= unique(arrayfun(@(x) x.tCond,cLabelOut,'UniformOutput',false));
stimGout= unique(arrayfun(@(x) x.stimGroup,cLabelOut,'UniformOutput',false));
%% update tempDat 
%1) fix names for tCond and stimGroup

    
    %Add numstim,pwr and tCond to the tempDat structure array.
    for i=1:length(tempDat)
        tempDat(i).tCond=cLabelOut(i).tCond;
        tempDat(i).stimGroup=cLabelOut(i).stimGroup;
    end
    
  
    % check my handi work
    tCond1= unique(arrayfun(@(x) x.tCond,tempDat,'UniformOutput',false));
%% upDate folder names from 2020-03-10 experiment.
% Here, all folder names should be 50ms exp rather than 5 and 50. This has
% been verified by the experimental script and the metaFrame.

dateInd=find(arrayfun(@(x) boolRegExp(x.date,'2020-03-10'), tempDat));

for i=1:length(dateInd)
    tempName=tempDat(dateInd(i)).fileName;
    pat='(?<=Exp_).*(?=_numStim)';
    repPat='50';
    IntPat='(?<=Int_).*(?=_Exp)';
    intVal=regexp(tempName,IntPat,'match');
    pwr=str2num(repPat)*str2num(intVal{1});
    newName=regexprep(tempName,pat,repPat);
    tempDat(dateInd(i)).fileName=newName;
    tempDat(dateInd(i)).dataSubdir=newName;
    tempDat(dateInd(i)).pwr=sprintf('%s',num2str(pwr));
    tempDat(dateInd(i)).stimGroup=sprintf('%s-%s-%s',tempDat(dateInd(i)).tCond,tempDat(dateInd(i)).pwr,tempDat(dateInd(i)).numStim);
    
end
%%
saveRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\Analyzed Data\2020-09-30_PLB-Cdc42-MultiStim';
fileName='2020-10-20-Cdc42-AllMultiPulse-v2.mat';

save([saveRoot filesep fileName],'tempDat','-v7.3');
%%
clear dat
%% load root
loadRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\Analyzed Data\2020-08-07-all-multistimDats';

    
%% load     
load([loadRoot filesep fileName]);

%% add stimGroups to dat and save
dat=tempDat;

%% save the dat
saveRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\2019-03-22_Cdc42Tk-Adaptaion-2Pulse exps\Analyzed Data\2019-04-03_PLB-PP1_Cdc42TK-2Pulse';
fileName='2021-03-28-Cdc42-All2Pulsedat-v2.mat';

save([saveRoot filesep fileName],'dat','-v7.3');
%% Make a data struct array (pDatPDL) that contains the relevant plot data
%2020-08-07: pDatPLD will be a struct array that has the length of the
%different conditons in tempDat.stimGroups. It will store the ratio and
%time data for each cells and the means. 

%parse exp condtions and preallocate pDatPDL
stimGroup= unique(arrayfun(@(x) x.stimGroup,dat,'UniformOutput',false));
stimGroup=natsortfiles(stimGroup);

%keepDate='2019-09-13';
%tCond=unique(arrayfun(@(x) x.tCond,dat,'UniformOutput',false));
for g=length(stimGroup):-1:1
    
    pDatPDL(g).stimGroup=[];
    pDatPDL(g).tCond=[];
    pDatPDL(g).pwr=[];
    pDatPDL(g).numStim=[];
    pDatPDL(g).sd=[];
    pDatPDL(g).sem=[];
    pDatPDL(g).meanNormRatio=[];
    pDatPDL(g).meanTKtime=[];
    pDatPDL(g).meanTktimeAdj=[];
    pDatPDL(g).meanStimTime=[];
    pDatPDL(g).meanStimTimeAdj=[];
    pDatPDL(g).t0isStimTime=[];
    pDatPDL(g).stimInStimRef=[];
    pDatPDL(g).allNormRatio=[];
    pDatPDL(g).tkTime=[];
    pDatPDL(g).stimTime=[];
    pDatPDL(g).fNames=[];
    pDatPDL(g).date=[];
    pDatPDL(g).numReplicates=[];
    
end


numBefore=10;
for stimG=1:length(stimGroup)
    clear allNormRatio tkTime stimTime fNames dates
    
    pDatPDL(stimG).stimGroup=stimGroup{stimG};
    %ind=find(arrayfun(@(x) boolRegExp(x.stimGroup,stimGroup{stimG}),dat));
    ind=find(arrayfun(@(x) strcmp(x.stimGroup,stimGroup{stimG}),dat));
    % break stimGroup into parts for easier name segmentation
    tempCell=regexp(stimGroup{stimG},'(.*\w\D\>)','match');
    pwr=regexp(stimGroup{stimG},'-(\d{0,5})-','match');
    numStim=regexp(stimGroup{stimG},'-(\d{0,2}$)','match');
    pDatPDL(stimG).tCond=tempCell;
    pDatPDL(stimG).pwr=pwr;
    pDatPDL(stimG).numStim=numStim;
    
    %remove specific dates from ind
%     keepInd=[];
%     for i=1:length(ind)
%        % dat(ind(i)).date
%        keepInd(i)=boolRegExp(dat(ind(i)).date,rmDate);
%     end
%     ind=ind(logical(keepInd));
    %preallocate allNormRatio & tkTime
    maxLpixMeanA=max(arrayfun(@(x) length(x.pix_meanA),dat));
    allNormRatio=nan(length(ind),maxLpixMeanA);
    
    maxLtkTime=max(arrayfun(@(x) length(x.time.tkTime),dat));
    tkTime=nan(length(ind),maxLtkTime);
    tkTimeAdj=nan(length(ind),maxLtkTime);
    
    maxLStimTime=max(arrayfun(@(x) length(x.time.frapTime),dat));
    stimTime=nan(length(ind),maxLStimTime);
    stimTimeAdj=nan(length(ind),maxLStimTime);
    
    for i=1:length(ind)
       % pDatPDL(stimG).datSubdir{i}=dat(ind(i)).fileName;
        
        tempPixMeanA=dat(ind(i)).pix_meanA./nanmean(dat(ind(i)).pix_meanA(:,1:numBefore));
        allNormRatio(i,1:length(tempPixMeanA))=tempPixMeanA; %normalized to the mean of the frames before stim (numBefore).
        % define number of well replicates
        numReps=length(~isnan(allNormRatio(:,2)));
        
        tempTKTime=dat(ind(i)).time.tkTime;
        tkTime(i,1:length(tempTKTime))=dat(ind(i)).time.tkTime;
        
        tempTKTimeAdj=dat(ind(i)).time.tkTimeAdj;
        tkTimeAdj(i,1:length(tempTKTime))=tempTKTimeAdj;
        
        tempStimT=dat(ind(i)).time.frapTime;
        tempStimTAdj=dat(ind(i)).time.frapTimeAdj;
        if isempty(tempStimT)
             stimTime(i,:)=nan;
        else
        stimTime(i,1:length(tempStimT))=tempStimT;
       stimTimeAdj(i,1:length(tempStimT))=tempStimTAdj;
        end
        fNames{i}=dat(ind(i)).fileName;
        dates{i}=dat(ind(i)).date;
    end
    sd=nanstd(allNormRatio);
    sem=nanstd(allNormRatio)./sqrt(size(allNormRatio,1));
    meanNormRatio=nanmean(allNormRatio);
    meanTKtime=nanmean(tkTime);
    meanStimTime=nanmean(stimTime);
    meanTKtimeAdj=mean(tkTimeAdj,'omitnan');
    meanStimTimeAdj=mean(stimTimeAdj,'omitnan');
    
    
    pDatPDL(stimG).sd=sd;
    pDatPDL(stimG).sem=sem;
    pDatPDL(stimG).meanNormRatio=meanNormRatio;
    pDatPDL(stimG).meanTKtime=round(meanTKtime,2);
    pDatPDL(stimG).meanStimTime=round(meanStimTime,2);
    pDatPDL(stimG).meanTKtimeAdj=round(meanTKtimeAdj,2);
    pDatPDL(stimG).meanStimTimeAdj=round(meanStimTimeAdj,2);
    
    pDatPDL(stimG).t0isStimTime= round(pDatPDL(stimG).meanTKtime-pDatPDL(stimG).meanStimTime(1),2);
    pDatPDL(stimG).stimInStimRef=round(pDatPDL(stimG).meanStimTime-pDatPDL(stimG).meanStimTime(1),2);
    
    pDatPDL(stimG).allNormRatio=allNormRatio;
    pDatPDL(stimG).tkTime=tkTime;
    pDatPDL(stimG).stimTime=stimTime;
    pDatPDL(stimG).fNames=fNames;
    pDatPDL(stimG).date=dates;
    pDatPDL(stimG).numReplicates=numReps;
    
    
end
%% save pDatPDL data.
%2020-08-07: save this data for use with plotting in this script or in the
%modeling scripts.
pDatPDL20203024LongV1=pDatPDL;
%tempSaveRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\Analyzed Data\2020-08-07-all-multistimDats';
save([saveRoot filesep 'pDatPDL20203024LongV1.mat'],'pDatPDL20203024LongV1','-v7.3');
%% check stim times


for f=1:length(pDatPDL)
scattInd(f)=any(~isnan(pDatPDL(f).meanStimTime));
end

testR=[1:length(pDatPDL)];
testR=testR(scattInd);

figure; hold on
for f=1:length(testR)
    yInd=linspace(1, .01*length(testR),length(testR));
    xInd=~isnan(pDatPDL(testR(f)).meanStimTime);
    x=pDatPDL(testR(f)).meanStimTime;
    x=x(xInd);
    scatter(x,repmat(yInd(f),1,length(x)))

end
%% set a temporary script root
scriptRoot='C:\Users\George\Documents\temp Data Folder\Scripts\2021-02-05-SC590-PDL';
%% save a copy of this script in the scriptRoot

Version='SC590-R01';
FileNameAndLocation=matlab.desktop.editor.getActiveFilename;
%get just the script name. want to save in scriptRoot location
bslashInds=regexp(FileNameAndLocation,'\\');
fName=FileNameAndLocation(bslashInds(end)+1:end);
fName(end-1:end)=[];
newbackup=sprintf('%s\\%s-backup%s.m',scriptRoot,fName,Version);
A = exist(newbackup,'file')
if (A~=0)
warning('Backup already exists for the current version')
else
    copyfile(FileNameAndLocation,newbackup);
end

%% Plot data from pDatPDL
%2020-08-07: This script allows the user to plot all conditions stored in
%pDatPDL or select conditions.

%Note: since we are now plotting on the x axis using the stored time data
%the experimetns that did not receive stimulation are much faster because
%the filter cube was not switching. For future experiments I need to set
%the display to false so that the imaging can go faster.
% This was updated 2020-10-13.

allCellTypes=unique(arrayfun(@(x) x.tCond,dat,'UniformOutput',false));

cellType={'all'}; %{'^Ctrl-Ret','C10-Ret'}
power={'500-'};
numStim={'all'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);

%inds2plot=fliplr(inds2plot);
 inds2plot=[1 4 inds2plot];

colors=parula(length(inds2plot)+2);
fig=figure; hold on;


for i=1:length(inds2plot)
    
    drawShadedErrorRegionRMNans(pDatPDL(inds2plot(i)).meanTKtime,...
        pDatPDL(inds2plot(i)).meanNormRatio,pDatPDL(inds2plot(i)).sem,...
        'color',colors(i,:),'edgealpha',.4,'facealpha',.4); % SEM
plot(pDatPDL(inds2plot(i)).meanTKtime,pDatPDL(inds2plot(i)).meanNormRatio,'Color',colors(i,:),'LineWidth',2);
legends{i}=sprintf('%s-%s reps',pDatPDL(inds2plot(i)).stimGroup, num2str(pDatPDL(inds2plot(i)).numReplicates));
end

hold on;

%timeBars = buildTimeScaleBarsForPDLPlots(pDatPDL,inds2plot, colors);

timeBars = buildTimeScaleBarsForPDLPlots(pDatPDL,inds2plot, colors,'overlay',0,'miny',0.98);%,...
   % 'imfieldname','t0isStimTime','stimfieldname','stimInStimRef');

%time bars are normally plotted short->long. need to reverse order
%indOrder=fliplr([1:length(timeBars)]);
indOrder=fliplr(indOrder);
maxNumStim=max(arrayfun(@(x) abs(str2num(x.numStim{1})), pDatPDL(inds2plot)));
if maxNumStim>1
    for t=1:length(timeBars)
        rectangle('Position',[timeBars(indOrder(t)).x,timeBars(indOrder(t)).y,timeBars(indOrder(t)).w,timeBars(indOrder(t)).h],'FaceColor',timeBars(t).colors,'EdgeColor','none');
    end
else
    t=1;
    rectangle('Position',[timeBars(t).x,.975,timeBars(t).w,1.07-.975],'FaceColor',[0.5 0.5 0.5], 'EdgeColor','none')
end
legend(legends,'Interpreter','none','Fontsize',6);
ylabel('Cdc42 Activity (Mean FRET Ratio)','FontSize',8)
xlabel('Time (s)','FontSize',8)
xlim([-15 190])
 ylim([.965 1.06])
%yticks([.99:0.01:1.045])
%title(sprintf('%s-PDL',cellType{1}),'Interpreter','none');
%% print number of replicates in a table
cellType={'all'}; %{'^Ctrl-Ret','C10-Ret'}
power={'-5000-'};
numStim={'all'}
inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);
conds=cell(length(inds2plot),1);
numRep=zeros(length(inds2plot),1);
clear tb;
for i=1:length(inds2plot)
    conds{i}=pDatPDL(inds2plot(i)).stimGroup;
    numRep(i)=pDatPDL(inds2plot(i)).numReplicates;
end

tb=table(conds,numRep)
%% plot all curves from a single condition in pDatPDL
ind=43
colors=parula(size(pDatPDL(ind).allNormRatio,1));
figure; hold on;


for i=1:size(pDatPDL(ind).allNormRatio,1)
    plot(pDatPDL(ind).meanTKtime,pDatPDL(ind).allNormRatio(i,:),'Color',colors(i,:),'LineWidth',2);
   
    fprintf('%s: %s\n',pDatPDL(ind).fNames{i},pDatPDL(ind).date{i}{1});
    
end
