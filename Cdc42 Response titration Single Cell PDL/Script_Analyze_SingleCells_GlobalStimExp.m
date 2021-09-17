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

%2020-12-18 Updated for sinle cell analysis
%2020-01-01 Updated to use  mean(fRatio(1:10)) <=deadThresh | isnan(mean(fRatio(1:10))) | any(yfp)>6e4 | any(cfp)>6e4
% as the new dead filter.
%% set roots

masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\';
% masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\GlobalStim-Long\';
%masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\Global 1 Stim All\';
%masterRoot = 'C:\Users\George\Documents\temp Data Folder\'; 
%% camera and donut corrections
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat',...
    'donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');

%ratioCorrectionImage
%load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\2020-11-12_20xTK RatioCorrection-post epi-Tirf mirror install\2D ratio correction-TK20xMaster.mat');

% load empty BG images
load(['D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\'...
    'global multiStim intensity titration\Analyzed Data\2020-02-04_PLB_Cdc42KO-PakKo-multistim_MedPwr'...
    filesep 'emptyP384WellsTK.mat'],'emptyYFP', 'emptyCFP');
 %% define the size of rawDataDir
% rawDataDir=getDirectories([masterRoot 'Raw Data'], 'SingleStim');
% %rawDataDir=rawDataDir(end);
 %% define cLabel, a structure array that contains label info for each cell
% % the struct will contain the cell treatment, number of stims and light
% % stim power data. The structure is defined outside of the primary data
% % processing loop so that all image folders can be identified across
% % multiple data directories. The total number of cells will also be used to
% % preallocate the dat variable, which is a struct array that contains the
% % processed image data.
% 
% 
% cLabel=[];
% for numRawDir=1:length(rawDataDir)
%     clear dataSubdir;
%     parentFolder=rawDataDir{numRawDir};
%     rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
%     dataSubdir=getSubdirectories(rawDataroot,'SC344'); % get the subdirectories for each folder in rawData
%     dataSubdir=natsortfiles(dataSubdir); % sort the data
%     
%     cellLine='PLB';
% %     constructLabels='SC344-SC147';
%     cTemp=getUniqueLabelsPDLAssay(dataSubdir,cellLine);
%     cLabel=[cLabel cTemp];
%     
% end
% 
 %% if label names are correct then make cLabelOut=cLabel
% cLabelOut=cLabel;
% %% adjust cLabels
% % multiStim label adj
% % sStr={'Cdc42KO','^C10$','^Ctrl$','^LatA$','^LatA-Ret$','Inh','LatA-C10','^LatA-Ctrl$'};
% % repStr={'C10-Ret','C10-Ret','Ctrl-Ret','LatA-Ctrl-Ret','LatA-Ctrl-Ret','PAK1Inhib-Ret','LatA-C10-Ret','LatA-Ctrl-Ret'};
% sStr={'10ugmlRetinal','^ret$','ret_LatA'};
% repStr={'Ctrl-Ret','Ctrl-Ret', 'LatA-Ctrl-Ret'};
% 
% cLabelOut=cLabelAdjPDLassay(cLabel,sStr,repStr);
% 
%     tCond1= unique(arrayfun(@(x) x.tCond,cLabelOut,'UniformOutput',false));
% stimGout= unique(arrayfun(@(x) x.stimGroup,cLabelOut,'UniformOutput',false));
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
    cTemp=getUniqueLabelsPDLAssay(dataSubdir,cellLine);
    cLabel=[cLabel cTemp];
    
end
%% adjust cLabels
sStr={'Cdc42KO','^C10$','^Ctrl$','^LatA$','^LatA-Ret$','Inh','LatA-C10','^LatA-Ctrl$'};
repStr={'C10-Ret','C10-Ret','Ctrl-Ret','LatA-Ctrl-Ret','LatA-Ctrl-Ret','PAK1Inhib-Ret','LatA-C10-Ret','LatA-Ctrl-Ret'};
cLabelOut=cLabelAdjPDLassay(cLabel,sStr,repStr);

    tCond1= unique(arrayfun(@(x) x.tCond,cLabelOut,'UniformOutput',false));
stimGout= unique(arrayfun(@(x) x.stimGroup,cLabelOut,'UniformOutput',false));

%% preallocate dat
for i=length(cLabel):-1:1
    dat(i).preMean1=[]
    dat(i).preMean2=[];
    dat(i).postMean=[];
    dat(i).singleFRET=[];
    dat(i).fileName=[];
    dat(i).singleYFP=[];
    dat(i).singleCFP=[];
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
    
    
    % select the correct ratio correction image
    %2D ratio correction-TK20xMaster.mat
    rawDataDate=regexp(rawDataDir{numRawDir},'[\d]+\-[\d]+\-[\d]+','match');
    date1=datetime(rawDataDate);
    dateThresh=datetime('2020-10-01');
    clear ratioCorrectionFinal;
    if date1>dateThresh
        load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\2020-11-12_20xTK RatioCorrection-post epi-Tirf mirror install\2D ratio correction-TK20xMaster.mat');
    else
        load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\2D ratio correction-TK20xMaster.mat');
    end
    
    
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
        filenames2=getFilenames(folder,imgName2) ;
        
        
        % get time stamps for the current set of fileNames
        if isempty(metaFrame)
            ts=buildTSforPDLdataPreMetaFiles(folder, filenames1);
        else
            ts = getImgTimesForCellMetaDataPDLassay(metaFrame, filenames1,numBefore);
        end
        
        filenames1=filenames1(1:17);
        filenames2=filenames2(1:17);
        
        cropsize=[951 921];
        imsY=cell(1,length(filenames1));
        imsC=cell(1,length(filenames2));
%         grayIm=cell(1,length(filenames1));
%         ratioImMovie=cell(1,length(filenames1));
        tic;
        for frameNum=1:length(filenames1)
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
%         density_scatter_heatmap(sumYDead+sumCDead,deadRatio,0:100:7e4,0:0.01:3);
%         hold on; yline(deadThresh,'-w')
        deadThresh=0.85;
        
        
        %preallocate for speed
        
        % whole cell preallocation:
        coors=cell(1,length(filenames2));
        pixList=cell(1,length(filenames2));
        mask=cell(1,length(filenames2));
        imYOut=cell(1,length(filenames2));
        imCout=cell(1,length(filenames2));
        for j=1:length(filenames2)
            tempCFP=imC(:,:,j); tempYFP=imY(:,:,j);
            [bgCFP,~]=imageSubtractBGScaleEmptyWell(tempCFP,pix,cropEyfp,cropEcfp,'cfp',true);
            [bgYFP,~]=imageSubtractBGScaleEmptyWell(tempYFP,pix,cropEyfp,cropEcfp,'yfp',true);
            imSum=bgCFP+bgYFP;
            [coors{j},pixList{j},mask{j}]=getNucleiPositions2(imSum,'minsize',100);
            
            %dont take the mean of the ratio image as small and large ratios
            %will bias the result. Instead take the mean of both channels
            %and then ratio.
            corrYFP=bgYFP./ratioCorrectionFinal;
            corrYFP(~mask{j})=nan;
            corrCFP=bgCFP;
            corrCFP(~mask{j})=nan;
            imYOut{j}=corrYFP;
            imCout{j}=corrCFP;
            
        end
        traj={}; rmInd=[]; rmInd=logical(rmInd);
        minDist=16;
        traj=basicTrackingFunction(coors,minDist);
        % filter traj to remove truncated data
        for t=1:length(traj)
            if numel(traj{t}) < length(filenames2)*9
                rmInd(t)=true;
            end
        end
        traj(rmInd)=[];
        dat(countSubdir).fileName=dataSubdir{subNum};
        tMask=cellfun(@(x) double(x), mask, 'UniformOutput', false);
        % for each traj, iterate through the pixel list
        for t1=1:length(traj) %iterate of traj
            %preallocate the fret vals below
            yfp=nan(1,length(filenames2));
            cfp=nan(1,length(filenames2));
            fRatio=nan(1,length(filenames2));
            for f=1:length(traj{t1}(:,9)) %count frames in image
                pixL=[];
                trajInd=traj{t1}(f,7);
                pixL=pixList{f}{trajInd};
                yfp(f)=nanmean(imYOut{f}(pixL));
                cfp(f)=nanmean(imCout{f}(pixL));
                fRatio(f)=yfp(f)/cfp(f);
                tMask{f}(pixL)=10; 
                
                
            end %f=1:length(traj{t1}(:,9))
            % deadThresh: set cells that have fret ratios <.9 to nan;
            if  mean(fRatio(1:10)) <=deadThresh | isnan(mean(fRatio(1:10))) | any(yfp)>6e4 | any(cfp)>6e4
                yfp=nan(1,length(filenames2));
                cfp=nan(1,length(filenames2));
                fRatio=nan(1,length(filenames2));
            end
                dat(countSubdir).singleCFP(t1,:)=cfp;
                dat(countSubdir).singleYFP(t1,:)=yfp;
                dat(countSubdir).singleFRET(t1,:)=fRatio;
                
                preMean1=nanmean(fRatio(1,1:5));
                preMean2=nanmean(fRatio(1,6:10));
                postMean=nanmean(fRatio(1,13:17));
                dat(countSubdir).preMean1(t1)=preMean1;
                dat(countSubdir).preMean2(t1)=preMean2;
                dat(countSubdir).postMean(t1)=postMean;
           
        end %for t1
       
        
        dat(countSubdir).fileName=dataSubdir{subNum};
        dat(countSubdir).date=regexp(rawDataDir{numRawDir}, '^[0-9\-]+', 'match');
        dat(countSubdir).time=ts;
        dat(countSubdir).tCond=cLabelOut(countSubdir).tCond;
        dat(countSubdir).numStim=cLabelOut(countSubdir).numStim;
        dat(countSubdir).pwr=cLabelOut(countSubdir).pwr;
        dat(countSubdir).dataSubdir=cLabelOut(countSubdir).dataSubdir;
        dat(countSubdir).stimGroup= sprintf('%s-%s-%s',cLabelOut(countSubdir).tCond,cLabelOut(countSubdir).pwr,cLabelOut(countSubdir).numStim);% =[tCond-pwr-numStim];
        
        
        %         moviePath=[saveRoot filesep dataSubdir{subNum} '.avi'];
        %         writeFRETMoviesFromImageSequences(ratioImMovie,grayIm,moviePath,...
        %             'imframerate',ts.frameRate,'objective',20,'boundsfret',[.8 1.2],'scalebarlength',75);
        toc;
        fprintf('\n');
        
    end % subNum=1:length(dataSubdir) ratioVals
    
end % numRawDir
%% save dat structure
SingleCellDat20201230MultiStim=dat;
save([saveRoot filesep 'SingleCellDat20201230MultiStim.mat'],'SingleCellDat20201230MultiStim','-v7.3');
%% save a copy of this script in the scriptRoot
scriptRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\Global 1 Stim All\Analyzed Data\2020-10-31-PLB-Cdc42-SingleStim';


Version='2';
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


%% analyze pseudo single cell data and create histogram data
novelStimG=unique(arrayfun(@(x) x.stimGroup, dat,'UniformOutput', false));

%preallocate
for i=length(novelStimG):-1:1
    hDat(i).allPeak=[];
    hDat(i).allPre=[];
    hDat(i).stimG=[];
    hDat(i).tCond=[];
    hDat(i).pwr=[];
    hDat(i).numStim=[];
    hDat(i).fretR=[];
end

for n=1:length(novelStimG)
    pat=sprintf('^%s',novelStimG{n});
    sgInds=find(arrayfun(@(x) boolRegExp(x.stimGroup, pat ), dat));
    hDat(n).stimG=novelStimG{n};
    for s=1:length(sgInds)
        if s==1
            hDat(n).tCond=dat(sgInds(s)).tCond;
            hDat(n).pwr=dat(sgInds(s)).pwr;
            hDat(n).numStim=dat(sgInds(s)).numStim;
        end
        ctrlRatio=dat(sgInds(s)).preMean1./dat(sgInds(s)).preMean2;
        peakRatio=dat(sgInds(s)).postMean./dat(sgInds(s)).preMean2;
        % remove nans
        crNanInds=isnan(ctrlRatio);
        prNanInds=isnan(peakRatio);
        
        ctrlRatio(crNanInds)=[];
        peakRatio(prNanInds)=[];
        
        allFret=[];
        allFret=dat(sgInds(s)).singleFRET;
        %rm nans for allFret
        afInds=[];
        for af=1:size(allFret,1)
            afInds(af)=all(isnan(allFret(af,:)));
        end 
        afInds=logical(afInds);
        allFret(afInds,:)=[];

        hDat(n).fretR=[hDat(n).fretR; allFret];
        hDat(n).allPre=[hDat(n).allPre ctrlRatio];
        hDat(n).allPeak=[hDat(n).allPeak peakRatio];
    end
end
%% filter hDat to see how different filters look
deadT=.925;
for h=1:length(hDat)
    filtTF=false(1,size(hDat(h).fretR,1));
    for f=1:size(hDat(h).fretR,1)
        filtTF(f)=mean(hDat(h).fretR(f,1:10))<deadT;
    end
    hDat(h).fretR(filtTF,:)=[];
    hDat(h).allPeak(filtTF)=[];
    hDat(h).allPre(filtTF)=[];
    
end

%%
fileName='2020-12-26-Cdc42-1Pulse-SingleCell-datAndhDat-v1.mat';

save([saveRoot filesep fileName],'dat','hDat','-v7.3');
%% plot histograms
allCellTypes=unique(arrayfun(@(x) x.tCond,hDat,'UniformOutput',false));

cellType={'^Ctrl-Ret'};
power={'50$','500$','5000$'};
numStim={'^0$','30'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(hDat,cellType,power,numStim);

color=parula(length(inds2plot)+1);
for i=1:length(color)
   c{i}=color(i,:); 
end 
%fig=figure; hold on;
legends={};
nbins=50;
fPanels=2
fNames=fieldnames(hDat);
fInds=boolRegExp(fNames,'all');
fNames=fNames(fInds);

for f=1:fPanels
    
    for h=1:length(inds2plot)
        hold on;
        s=subplot(1,2,f);
       %histogram(vect(hDat(inds2plot(h)).(fNames{f})),nbins,'FaceAlpha',0.5,'EdgeColor',color(h,:));
                %histogram(vect(hDat(inds2plot(h)).(fNames{f})),nbins,'smooth');
      plotFrequencyDistributions(hDat(inds2plot(h)).(fNames{f}), 0.95:0.005:1.15,c(h),0);
      
      
      legends{h}=sprintf('%s-n=%s',hDat(inds2plot(h)).stimG,num2str(length(hDat(inds2plot(h)).(fNames{f}))));
      
    end
    title(sprintf('%s',fNames{f}));
    ylabel('Fraction of Cells');
    xlabel('FRET Ratio Fold Change');
end
legend(legends);

%xlim([0.8 1.3])
%title('Cdc42-KO');
%% plot all fret ratios from a trajectroy
allCellTypes=unique(arrayfun(@(x) x.tCond,hDat,'UniformOutput',false));

cellType={'^C10-Ret'};
% power={'^0','^1000$','^5000$'};
% numStim={'1'}; %{'-0$','-1$','-7','-15','-30'}
power={'50$','500$','5000$'};
numStim={'^0$','30'}; %{'-0$','-1$','-7','-15','-30'}


inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(hDat,cellType,power,numStim);
xVal=[1:17];
for s=1:length(inds2plot)
    subplot(1,length(inds2plot),s); hold on;
    for h=1:size(hDat(inds2plot(s)).fretR,1)
        if h<1864
         plot(xVal,hDat(inds2plot(s)).fretR(h,:)./nanmean(hDat(inds2plot(s)).fretR(h,1:10)));
        end
    end
    ylabel('Norm. Single Cell Mean Fret Ratio');
    xlabel('Frame');
    title(sprintf('%s',hDat(inds2plot(s)).stimG));
end

%% violin plots
allCellTypes=unique(arrayfun(@(x) x.tCond,hDat,'UniformOutput',false));

cellType={'^Ctrl'};
power={'^0','^1000$','^5000$'};
pLabel={'0','1000','5000'};
numStim={'1'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(hDat,cellType,power,numStim);

color=parula(length(power)+1);
for i=1:length(color)
   c{i}=color(i,:); 
end 
%fig=figure; hold on;
legends={};
nbins=50;
fPanels=2
fNames=fieldnames(hDat);
fInds=boolRegExp(fNames,'all');
fNames=fNames(fInds);

violinCell=cell(1,length(inds2plot));

for f=1%:fPanels
    
    for h=1:length(inds2plot)
        hold on;
       violinCell{h}=hDat(inds2plot(h)).(fNames{f})';
    end
    %violin(violinCell);
    violin(violinCell,'xlabel',pLabel);
    title(sprintf('%s',cellType{1}));
    ylabel('Mean Fret Ratio')
    xlabel('Stimulation Power')
end
legend(legends)

%% add scatter to violin plot

hold on;
yd=[];
xd=[];

for f=1:2
    color={'b','r'};
   for h=1:length(inds2plot)
        
       yd=hDat(inds2plot(h)).(fNames{f})';
       xd=(ones(size(yd))*h).*(1+(rand(size(yd))-0.75)/10);
       %scatter(xd,yd,color{f},'filled', 'MarkerFaceAlpha',0.4);
       scatter(xd,yd,color{f},'filled', 'MarkerFaceAlpha',0.4);


    end
end
%xlim([0.8 1.3])
%title('Cdc42-KO
%% scatter allPeak v allPre

cellType={'^C10-Ret'};
power={'^5000$'};
numStim={'30'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(hDat,cellType,power,numStim);

hold on;
allPre1=[];
allPeak1=[];

color='b'
for h=1:length(inds2plot)
    
    allPre1=hDat(inds2plot(h)).allPre';
    allPeak1=hDat(inds2plot(h)).allPeak';
   
    scatter(allPre1,allPeak1,color,'filled', 'MarkerFaceAlpha',0.4);
    
end

%% plot CDFs
% use this to estimate the population of responders v nonResponders

cellType={'^C10'};
power={'^0','^1000$','^5000$'};
numStim={'1'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(hDat,cellType,power,numStim);
nbins=50;
fPanels=2
fNames=fieldnames(hDat);
fInds=boolRegExp(fNames,'all');
fNames=fNames(fInds);



for h=1:length(inds2plot)
    
    allPre1=[];
    allPeak1=[];
    
    allPre1=hDat(inds2plot(h)).allPre';
    allPeak1=hDat(inds2plot(h)).allPeak';
    
    plotCDFs([allPre1,allPeak1]);
    
    legends{h}=sprintf('%s-n=%s',hDat(inds2plot(h)).stimG,num2str(length(hDat(inds2plot(h)).(fNames{f}))));
    
    
    title(sprintf('%s',hDat(inds2plot(h)).stimG));
    ylabel('Fraction of Cells')
    xlabel('FRET Ratio Fold Change')
end

%% density scatter meanPreStim v allPeak

cellType={'^C10-Ret'};
power={'50$','500$','5000$'};
numStim={,'^0','30'}; %{'-0$','-1$','-7','-15','-30'}

inds2plot=[];
inds2plot=sortpDatPLBbyStimConds(hDat,cellType,power,numStim);



for h=1:length(inds2plot)
    sPlot=subplot(2,2,h); hold on;
    allPreMean=[];
    for i=1:size(hDat(inds2plot(h)).fretR,1)
        allPreMean(i)=nanmean(hDat(inds2plot(h)).fretR(i,1:10));
    end
    allPeak1=[]
    allPeak1=hDat(inds2plot(h)).allPeak';
    density_scatter_heatmap(allPreMean,allPeak1,0.9:0.005:1.2,0.9:0.005:1.2,sPlot);
    title(sprintf('%s',hDat(inds2plot(h)).stimG));
    ylabel('allPeak1');
    xlabel('allPreMean');
end


