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

%2021-03-29 Modified to make movies for the PDL assays. I commented out
%most of the single cell processing steps. Additionally, I altered the
%inputs to getNucleiPositions2 so that imclearborder was turned off, and
%the min size was 50, max size 15000. This reduces the number of cells
%entering and leavign the masks between frames.

%2021-04-08 GB Updated to make movies that are multi paneled. I tailored
%this script for finding 3 specific files from one folder. I may need to
%change that to make the 1 pulse pdl figs. 
%% set roots
%masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\global multiStim intensity titration\';
% masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\GlobalStim-Long\';
masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\Global 1 Stim All\';
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
%rawDataDir=getDirectories([masterRoot 'Raw Data'], '2020-03-03_PLB_Cdc42KO-PakKo-multistim_MedPwr');
rawDataDir=getDirectories([masterRoot 'Raw Data'], '2018-05-17_PP1-Cdc42TK-TitrateStimInt-PLB-SingleStim');
%rawDataDir=getDirectories([masterRoot 'Raw Data'], '2020-03-05_PLB_C10_multistim');

%rawDataDir=rawDataDir(end);
 % define cLabel, a structure array that contains label info for each cell
% the struct will contain the cell treatment, number of stims and light
% stim power data. The structure is defined outside of the primary data
% processing loop so that all image folders can be identified across
% multiple data directories. The total number of cells will also be used to
% preallocate the dat variable, which is a struct array that contains the
% processed image data.
% cellFind={'A12_PLB-SC344-SC147-Ctrl-NoR_Int_10_Exp_50_numStims_30','B07_PLB-SC344-SC147-Ctrl-Ret_Int_10_Exp_50_numStims_0',...
%     'B08_PLB-SC344-SC147-Ctrl-Ret_Int_10_Exp_50_numStims_30'}; %for sMovie 1;

cellFind={'B02_PLB_SC344_SC147_10ugmlRetinal_Int_0_Exp_50','B03_PLB_SC344_SC147_10ugmlRetinal_Int_20_Exp_50',...
    'C04_PLB_SC344_SC147_10ugmlRetinal_Int_100_Exp_50'};

cLabel=[];
for numRawDir=1:length(rawDataDir)
    clear dataSubdir;
    parentFolder=rawDataDir{numRawDir};
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    dataSubdir=getSubdirectories(rawDataroot,'SC344'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the data
  
    for c=1:length(cellFind)
    ind(c)=find(boolRegExp(dataSubdir,cellFind{c}));
    end    
    
    dataSubdir=dataSubdir(ind);
    cellLine='PLB';
%     constructLabels='SC344-SC147';
    cTemp=getUniqueLabelsPDLAssay(dataSubdir,cellLine);
    cLabel=[cLabel cTemp];
    
end

 % if label names are correct then make cLabelOut=cLabel
cLabelOut=cLabel;
%% adjust cLabels
% multiStim label adj
% sStr={'Cdc42KO','^C10$','^Ctrl$','^LatA$','^LatA-Ret$','Inh','LatA-C10','^LatA-Ctrl$'};
% repStr={'C10-Ret','C10-Ret','Ctrl-Ret','LatA-Ctrl-Ret','LatA-Ctrl-Ret','PAK1Inhib-Ret','LatA-C10-Ret','LatA-Ctrl-Ret'};
sStr={'10ugmlRetinal','^ret$','ret_LatA'};
repStr={'Ctrl-Ret','Ctrl-Ret', 'LatA-Ctrl-Ret'};

cLabelOut=cLabelAdjPDLassay(cLabel,sStr,repStr);

    tCond1= unique(arrayfun(@(x) x.tCond,cLabelOut,'UniformOutput',false));
stimGout= unique(arrayfun(@(x) x.stimGroup,cLabelOut,'UniformOutput',false));


%% preallocate dat
for i=length(cLabel):-1:1
    mDat(i).imGray=[];
    mDat(i).imFRET=[];
    mDat(i).fileName=[];
    mDat(i).date=[];
    mDat(i).time=[];
    mDat(i).tCond=[];
    mDat(i).numStim=[];
    mDat(i).pwr=[];
    mDat(i).dataSubdir=[];
    mDat(i).stimGroup=[]; % =[tCond-pwr-numStim];
end
%% compute alignment and ratio correction data for each experimental day


countSubdir=0; % counter for total number of subdirectories read, allows me to build a big struct for 3 days of exp
TFtomKat=true;
 
imgName1='TomKat-1';
imgName2='TomKat-2'; 
%cropsize=[951 921];
cropsize=[750 750];
for numRawDir=1:length(rawDataDir)
    clear dataSubdir;
    parentFolder=rawDataDir{numRawDir};
    %     mkdir([masterRoot 'Scripts' filesep parentFolder]); % makes a folder for the exp day in the scripts folder
    %     mkdir([masterRoot 'Analyzed Data' filesep parentFolder]);% makes a folder for the exp day in the Analyzed Data folder
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    dataSubdir=getSubdirectories(rawDataroot,'SC344'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the data
    
    %pick specific folders
    for c=1:length(cellFind)
    ind(c)=find(boolRegExp(dataSubdir,cellFind{c}));
    end    
    dataSubdir=dataSubdir(ind);
    
    
    
    
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
    
    ratioCorrectionFinal=cropImMidOut(ratioCorrectionFinal,'cropsize', cropsize);
    
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
        
%         filenames1=filenames1(1:17);
%         filenames2=filenames2(1:17);
        
        %cropsize=[951 921];
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
        ratioImMovie=cell(1,length(filenames2));
        grayIm=cell(1,length(filenames2));
        
        for j=1:length(filenames2)
            tempCFP=imC(:,:,j); tempYFP=imY(:,:,j);
            [bgCFP,~]=imageSubtractBGScaleEmptyWell(tempCFP,pix,cropEyfp,cropEcfp,'cfp',true);
            [bgYFP,~]=imageSubtractBGScaleEmptyWell(tempYFP,pix,cropEyfp,cropEcfp,'yfp',true);
            imSum=bgCFP+bgYFP;
            %[coors{j},pixList{j},mask{j}]=getNucleiPositions2(imSum,'minsize',100);%
            %for single cells
            
            %movie version. less strict size requirements so that cells are
            %not coming into and out of the mask
            [coors{j},pixList{j},mask{j}]=getNucleiPositions2(imSum,'minsize',50,'maxsize',15000,...
                'imclearborder',false);
            
            %dont take the mean of the ratio image as small and large ratios
            %will bias the result. Instead take the mean of both channels
            %and then ratio.
            corrYFP=bgYFP./ratioCorrectionFinal;
            corrYFP(~mask{j})=nan;
            corrCFP=bgCFP;
            corrCFP(~mask{j})=nan;
            imYOut{j}=corrYFP;
            imCout{j}=corrCFP;
            ratioImMovie{j}=corrYFP./corrCFP;
            grayIm{j}=imSum;
            
        end
        %
        mDat(countSubdir).time=ts;
        mDat(countSubdir).imGray=grayIm;
        mDat(countSubdir).imFRET=ratioImMovie;
        mDat(countSubdir).fileName=folder;
        mDat(countSubdir).date=regexp(rawDataDir{numRawDir}, '^[0-9\-]+', 'match');
        mDat(countSubdir).tCond=cLabelOut(countSubdir).tCond;
        mDat(countSubdir).numStim=cLabelOut(countSubdir).numStim;
        mDat(countSubdir).pwr=cLabelOut(countSubdir).pwr;
        mDat(countSubdir).dataSubdir=cLabelOut(countSubdir).dataSubdir;
        mDat(countSubdir).stimGroup=sprintf('%s-%s-%s',cLabelOut(countSubdir).tCond,cLabelOut(countSubdir).pwr,cLabelOut(countSubdir).numStim);% =[tCond-pwr-numStim];
        
        toc;
        fprintf('\n');
        
    end % subNum=1:length(dataSubdir) ratioVals
    
end % numRawDir
%%
mDatSave=mDat;
%% save dat structure
SingleCellDat20210107SinlgeStim=mDat;
save([saveRoot filesep 'SingleCellDat20210107SinlgeStim.mat'],'SingleCellDat20210107SinlgeStim','-v7.3');
%% build a multipanel movie from mDat

 %make dat with ts. just to add the correct time stamps
 for m=1:length(mDat)
        mDat(m)=adjTimeStampsWithEmpiricalTimeMeasurements(mDat(m));
%         temptkTime=nan(1,length(mDat(countSubdir).time.tkTimeAdj));
%         temptkTime=round(mDat(countSubdir).time.tkTimeAdj,1);
%         
%         tempStimTime=nan(1,length(mDat(countSubdir).time.frapTimeAdj));
%         tempStimTime=mDat(countSubdir).time.frapTimeAdj;
 end
        moviePath=[saveRoot filesep 'sMovie3-multipanel.avi'];
%         writeFRETMoviesFromImageSequences(ratioImMovie,grayIm,moviePath,...
%             'tstktime',temptkTime,'objective',20,'boundsfret',[.99 1.1],'scalebarlength',75,...
%             'scalebartexttf',false,'tsstimtime',tempStimTime);

%extract data from dat to input into the function;
% imDim=[size(mDat(1).imFRET,1), size(mDat(1).imFRET,2), numel(mDat)];
% 
% imGrayAll=cell(imDim);
% imRatioAll=cell(imDim); 
labels={'No UV', 'PWR=20', 'PWR=100'};
writeFRETMoviesFromImageSequencesMultiPanel(mDat,moviePath,...
           'objective',20,'boundsfret',[.99 1.1],'scalebarlength',75,...
            'scalebartexttf',false,'labels',labels);
%% save a copy of this script in the scriptRoot

Version='sMovie1';
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


