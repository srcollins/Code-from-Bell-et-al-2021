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

%2021-03-28 Update: Updated the pDatPDL block to include more accurate time
%stamps for premetadata images.

%This script is specific for plotting the 2-pulse stimulation assays from
%Fig 3a&b and SFig2
%% set roots

masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\2019-03-22_Cdc42Tk-Adaptaion-2Pulse exps\';
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
rawDataDir=getDirectories([masterRoot 'Raw Data'], 'PLB');
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
    dataSubdir=getSubdirectories(rawDataroot,'SC344'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the data
    
    cellLine='PLB';
%     constructLabels='SC344-SC147';
    cTemp=getUniqueLabelsPDLAssay(dataSubdir,cellLine,'numstim','(?<=secPostStim1_)[0-9]{0,3}');
    cLabel=[cLabel cTemp];
    
end
%% adjust cLabels
sStr={''};
repStr={'Ctrl-Ret'};
cLabelOut=cLabelAdjPDLassay(cLabel,sStr,repStr,'all',true);


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
fileName='2020-11-14-Cdc42-All2Pulse-v1.mat';

save([saveRoot filesep fileName],'dat','-v7.3');


%% store dat as a new variable so that I can load the old dat to 
%concatenate the new data into the total dat structure.
oldestDat=dat;
clear dat

%% adjust cLabels
% block standardizes saved stimulation condition names for downstream use.
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
%% load root
loadRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\Analyzed Data\2020-08-07-all-multistimDats';

    
%% load     
load([loadRoot filesep fileName]);

%% add stimGroups to dat and save
dat=tempDat;

%% add adjusted time stamps for 2Pulse exp so that the time stamps onthe plots are correct
dat=adjTimeStampsWithEmpiricalTimeMeasurements(dat);
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
pDatPDL2Pulse20210107v1=pDatPDL;
tempSaveRoot='C:\Users\George\Dropbox\GB-PP1 manuscript\Opsin Paper figure Data\2-Pulse Fig Data\Processed Data';
save([tempSaveRoot filesep 'pDatPDL2Pulse20210107v1.mat'],'pDatPDL2Pulse20210107v1','-v7.3');
%% load
load([saveRoot filesep 'pDatPDL20201022v1.mat']);
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

%% Plot data from pDatPDL
%2020-08-07: This script allows the user to plot all conditions stored in
%pDatPDL or select conditions.

%Note: since we are now plotting on the x axis using the stored time data
%the experimetns that did not receive stimulation are much faster because
%the filter cube was not switching. For future experiments I need to set
%the display to false so that the imaging can go faster.
% This was updated 2020-10-13.

allCellTypes=unique(arrayfun(@(x) x.tCond{1},pDatPDL,'UniformOutput',false));

cellType={'all'};
power={'-500-'};
numStim={'-5','26','-50'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);

colors=parula(length(inds2plot)+1);
fig=figure; hold on;


for i=1:length(inds2plot)
    subplot(3,1,i); hold on;
    drawShadedErrorRegionRMNans(pDatPDL(inds2plot(i)).meanTKtime,pDatPDL(inds2plot(i)).meanNormRatio,pDatPDL(inds2plot(i)).sem,'color',colors(i,:)); % SEM
plot(pDatPDL(inds2plot(i)).meanTKtime,pDatPDL(inds2plot(i)).meanNormRatio,'Color',colors(i,:),'LineWidth',2);
legends{i}=sprintf('%s-%s reps',pDatPDL(inds2plot(i)).stimGroup, num2str(pDatPDL(inds2plot(i)).numReplicates));
end

hold on;
timeBars = buildTimeScaleBarsForPDLPlots(pDatPDL,inds2plot, colors);
for t=1:length(timeBars)
 rectangle('Position',[timeBars(t).x,timeBars(t).y,timeBars(t).w,timeBars(t).h],'FaceColor',timeBars(t).colors,'EdgeColor','none');
end
legend(legends,'Interpreter','none');
ylabel('Mean FRET Ratio')
xlabel('time (sec)')
xlim([-15 100])
ylim([.975 1.07])
title(sprintf('%s-PDL',cellType{1}),'Interpreter','none');

%% plot all curves from a single condition in pDatPDL
ind=8
colors=parula(size(pDatPDL(ind).allNormRatio,1));
figure; hold on;
   
  
for i=1:size(pDatPDL(ind).allNormRatio,1)
    plot(pDatPDL(ind).meanTKtime,pDatPDL(ind).allNormRatio(i,:),'Color',colors(i,:),'LineWidth',2);
   
    fprintf('%s: %s\n',pDatPDL(ind).fNames{i},pDatPDL(ind).date{i}{1});
   
end
title(pDatPDL(ind).stimGroup)
ylabel('Mean Fret Ratio')
xlabel('Time (sec)')

%% collect max vals for each second stim
% the goal here is to pull out the second stim peak value and the time
% between stim value. I will use the second stim time to crop the data set
% so that the find peaks function can find the correct peak.
clear legends
%restrict the criteria: only max pwr stims
cellType={'all'};
power={'-500-'};
numStim={'all'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);
colors=parula(length(inds2plot)+1);

% Ive added tkTimeAdj with more accurate time stamps. Use the fieldname
% here to specify
timeFname='meanTKtimeAdj';
stimFname='meanStimTimeAdj';

for pd=length(inds2plot):-1:1
peakDat(pd).peakMax=[];
peakDat(pd).tBetween=[];
end

for i=1:length(inds2plot)
    tInd=min(find(pDatPDL(inds2plot(i)).(timeFname)>pDatPDL(inds2plot(i)).(stimFname)(2)));
    peakDat(i).peakMax=max(findpeaks(pDatPDL(inds2plot(i)).meanNormRatio(tInd:end)));
    peakDat(i).tBetween=round(pDatPDL(inds2plot(i)).(stimFname)(2)-pDatPDL(inds2plot(i)).(stimFname)(1));
end



% pull out the data for the specific 
cellType={'all'};
power={'-500-'};
numStim={'-30$'};
exInd2Plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);

figure; hold on;
exColor=[0 0 0]
drawShadedErrorRegionRMNans(pDatPDL(exInd2Plot).(timeFname),pDatPDL(exInd2Plot).meanNormRatio,pDatPDL(exInd2Plot).sem,'color',exColor); % SEM
plot(pDatPDL(exInd2Plot).(timeFname),pDatPDL(exInd2Plot).meanNormRatio,'Color',exColor,'LineWidth',2);
curveLeg{1}=sprintf('%s-%s reps',pDatPDL(exInd2Plot).stimGroup, num2str(pDatPDL(exInd2Plot).numReplicates));
hold on;
for j=1:length(peakDat)
   scatter(peakDat(j).tBetween, peakDat(j).peakMax,'MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:)); 
   legends{j}=sprintf('%s-tBetween', num2str(peakDat(j).tBetween));
end

allLegends=[curveLeg legends];
legend(allLegends)
ylabel('Mean FRET Ratio')
xlabel('Time (s)')

%% ratio peak2/peak1
%% collect max vals for each second stim
% the goal here is to pull out the second stim peak value and the time
% between stim value. I will use the second stim time to crop the data set
% so that the find peaks function can find the correct peak.

% used to plot figure 2f 20201217 GB;
% updated 2020-12-30 to collect the peaks from the population rather than
% the final mean. THis will allow me to add errorbars

clear legends
%restrict the criteria: only max pwr stims
cellType={'all'};
power={'-500-'};
numStim={'all'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);
%inds2plot(1)=[];
colors=parula(length(inds2plot)+1);

for pd=length(inds2plot):-1:1
peakDat(pd).peakMax1=[];  
peakDat(pd).peakMax2=[];
peakDat(pd).tBetween=[];
peakDat(pd).semFChange=[];
peakDat(pd).sdFChange=[];
peakDat(pd).fChange=[];
peakDat(pd).numRep=[];
end

peaks=zeros(1,2);
for i=1:length(inds2plot)
    fChange=zeros(1,size(pDatPDL(inds2plot(i)).allNormRatio,1));
    peakDat(i).numRep=size(pDatPDL(inds2plot(i)).allNormRatio,1);
    for subPeak=1:size(pDatPDL(inds2plot(i)).allNormRatio,1)
        peaks(subPeak,:)=findpeaks(pDatPDL(inds2plot(i)).allNormRatio(subPeak,:),'MinPeakHeight',1.01);
        peakDat(i).peakMax1(subPeak)=peaks(subPeak,1);
        
        if numel(peaks(subPeak,:))>1
            peakDat(i).peakMax2(subPeak)=peaks(subPeak,2);
        end
        fChange(subPeak)=(peakDat(i).peakMax2(subPeak)-1)/(peakDat(i).peakMax1(subPeak)-1);
        peakDat(i).tBetween=round(pDatPDL(inds2plot(i)).(stimFname)(2)-pDatPDL(inds2plot(i)).(stimFname)(1));
    end
    
    peakDat(i).fChange=mean(fChange);
    peakDat(i).semFChange=std(fChange)/sqrt(length(fChange));
    peakDat(i).sdFChange=std(fChange);
    
end



% pull out the data for the specific 
cellType={'all'};
power={'-500-'};
numStim={'-30$'};
exInd2Plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);

figure; hold on;
exColor=[0 0 0]
%drawShadedErrorRegionRMNans(pDatPDL(exInd2Plot).(timeFname),pDatPDL(exInd2Plot).meanNormRatio,pDatPDL(exInd2Plot).sem,'color',exColor); % SEM
%plot(pDatPDL(exInd2Plot).(timeFname),pDatPDL(exInd2Plot).meanNormRatio,'Color',exColor,'LineWidth',2);
%curveLeg{1}=sprintf('%s-%s reps',pDatPDL(exInd2Plot).stimGroup, num2str(pDatPDL(exInd2Plot).numReplicates));
hold on;

for j=1:length(peakDat)
    
  fig=errorbar(peakDat(j).tBetween,peakDat(j).fChange,peakDat(j).semFChange,'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
      %fig=scatter(peakDat(j).tBetween, (peakDat(j).peakMax2-1)/(peakDat(j).peakMax1-1),'MarkerFaceColor','k','MarkerEdgeColor','k'); 
fig.Color='k';
    
   legends{j}=sprintf('%s sec, n=%s', num2str(peakDat(j).tBetween),num2str(peakDat(j).numRep));
end

%allLegends=[curveLeg legends];
%legend(allLegends)
%legend(legends)

ylabel('Relative Amplitude of 2nd Peak')
xlabel('Time Between Stimuli (s)')
xlim([10 55]);
%%

%% Subplot 4 curves 
%2020-08-07: This script allows the user to plot all conditions stored in
%pDatPDL or select conditions.

%Note: since we are now plotting on the x axis using the stored time data
%the experimetns that did not receive stimulation are much faster because
%the filter cube was not switching. For future experiments I need to set
%the display to false so that the imaging can go faster.
% This was updated 2020-10-13.

%2020-12-15GB I used this block to create the 3 and 4 panel subfigures for
%the 2 pulse data. 4 panel fig was used for Fig 2e 20201217
legends={};
allCellTypes=unique(arrayfun(@(x) x.tCond,dat,'UniformOutput',false));

cellType={'all'};
power={'-500-'};
numStim={'-5','15','26','-50'}; 

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);

colors=parula(length(inds2plot)+1);
fig=figure; hold on;


for i=1:length(inds2plot)
    subplot(2,2,i); hold on;
    
    hold on;
minY=.975; % use the bottom of the ylimit for the axis
yaxismax=1.06; % use the max ylim for the axis
g=[0.5 0.5 0.5];
tBarColor=repmat(g,size(colors,1),1);
timeBars = buildTimeScaleBarsFor2PulsePDLPlots(pDatPDL,minY,inds2plot, tBarColor,'yaxismax',yaxismax,'stimfieldname','meanStimTimeAdj','tktimefieldname','meanTKtimeAdj');

    for j=1:2 %pulse number
 rectangle('Position',[timeBars(i).x(j),timeBars(i).y(j),timeBars(i).w(j),timeBars(i).h(j)],'FaceColor',timeBars(i).colors,'EdgeColor','none');
    end

    
    drawShadedErrorRegionRMNans(pDatPDL(inds2plot(i)).(timeFname),pDatPDL(inds2plot(i)).meanNormRatio,pDatPDL(inds2plot(i)).sem,'color',colors(i,:),'edgealpha',0); % SEM
    plot(pDatPDL(inds2plot(i)).(timeFname),pDatPDL(inds2plot(i)).meanNormRatio,'Color',colors(i,:),'LineWidth',1);
    %errorbar(pDatPDL(inds2plot(i)).(timeFname),pDatPDL(inds2plot(i)).meanNormRatio,pDatPDL(inds2plot(i)).sem,'Color',colors(i,:),'LineWidth',1.5);
    
   % tLabel=sprintf('Delay %s s', num2str(round(pDatPDL(subPGroups{s}(1)).meanStimTimeAdj(2),0)));

    tVals=round(abs(diff(pDatPDL(inds2plot(i)).meanStimTimeAdj(1,:))));
    title(sprintf('Delay %s s',num2str(tVals)));
    
    xlim([-15 100])
    ylim([.975 1.06])
    
    
    
end

ylabel('Cdc42 Activity (Mean FRET Ratio)')
    xlabel('Time (s)')
%legend(legends,'Interpreter','none');

%title(sprintf('%s-PDL',cellType{1}),'Interpreter','none');
%% Subplot 2 pulse PDL exp data
%Used to make SFig2
allCellTypes=unique(arrayfun(@(x) x.tCond,dat,'UniformOutput',false));
allStimDelays=arrayfun(@(x) regexp(x.fileName,'(?<=secPostStim1).[\d]+','match'),dat,'UniformOutput',false);
stimDelays=unique(cellfun(@(x) x,allStimDelays));
stimDelays=natsortfiles(stimDelays);
cellType={'^Ctrl'}; %{'^Ctrl-Ret','C10-Ret'}
power={'125','250','375','500'};
numStim={'all'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(pDatPDL,cellType,power,numStim);

subPGroups=cell(1,length(stimDelays));
% sort by secPostStim1
for i=1:length(inds2plot)
    tempInd=inds2plot(i);
    tempfName=regexp(pDatPDL(inds2plot(i)).fNames{1},'(?<=secPostStim1).[\d]+','match');
    tempfName=sprintf('%s$',tempfName{1});
    groupInd=find(boolRegExp(stimDelays,tempfName));
    subPGroups{groupInd}=[subPGroups{groupInd},tempInd];
end
%inds2plot=fliplr(inds2plot);
%inds2plot=[24 inds2plot];
fg=figure('Units','inches','Position',[-.25,-.25,8.5,11]);
width=3;


colors=parula(4);
for s=1:length(subPGroups)
    hold on;
    subplot(ceil(length(subPGroups)/width),width,s);
  
   
    
     timeBars=buildTimeScaleBarsFor2PulsePDLPlots(pDatPDL,0.985,subPGroups{s}, colors,'stimfieldname','meanStimTimeAdj','tktimefieldname','meanTKtimeAdj');
    %time bars are normally plotted short->long. need to reverse order
   
    for t=1:length(timeBars)
        
        rectangle('Position',[timeBars(t).x(1),timeBars(t).y(1),timeBars(t).w(1),timeBars(t).h(1)],'FaceColor',[0.5 0.5 0.5], 'EdgeColor','none');
        rectangle('Position',[timeBars(t).x(2),timeBars(t).y(2),timeBars(t).w(2),timeBars(t).h(2)],'FaceColor',[0.5 0.5 0.5], 'EdgeColor','none');
    end
    
    legends={1,length(subPGroups{s})};
    for i=1:length(subPGroups{s})
        hold on;
        drawShadedErrorRegionRMNans(pDatPDL(subPGroups{s}(i)).meanTKtimeAdj,...
            pDatPDL(subPGroups{s}(i)).meanNormRatio,pDatPDL(subPGroups{s}(i)).sem,...
            'color',colors(i,:),'edgealpha',0,'facealpha',.4); % SEM
        plot(pDatPDL(subPGroups{s}(i)).meanTKtimeAdj,pDatPDL(subPGroups{s}(i)).meanNormRatio,'Color',colors(i,:),'LineWidth',1);
        pwrTemp={}; pwrTemp=regexp(pDatPDL(subPGroups{s}(i)).pwr,'[\d]+','match');
        pwr=str2double(pwrTemp{1})/50;
        legends{i}=sprintf('%s',num2str(pwr));
        %legends{i}=sprintf('%s-%s reps',pDatPDL(subPGroups{s}(i)).stimGroup, num2str(pDatPDL(subPGroups{s}(i)).numReplicates));
    end
    
    hold on;
   
    tLabel=sprintf('Delay %s s', num2str(round(pDatPDL(subPGroups{s}(1)).meanStimTimeAdj(2),0)));
    title(tLabel,'FontSize',8)
    legend(legends,'Interpreter','none','FontSize',6);
    ylabel('Cdc42 Activity (Mean FRET Ratio)','FontSize',8)
    xlabel('Time (s)','FontSize',8 )
    xlim([-15 120])
    ylim([.985 1.055])
end % end for s
