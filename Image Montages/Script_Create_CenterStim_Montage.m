%% Script to compute cell montages

%THis script was used to make the Centerstim montage (Fig 7a). 
%See bottom of script for making montage from loaded imdat. 

%This script will save processed image data for each cell in a structure called imDat.
%After all experiments are processed, each imDat is accessed individually
%for generating the montages.
%% Root for saving data
%masterRoot = 'C:\Users\George\Documents\temp Data Folder\';  % This should be the path to the folder containing all the data from the experiment
masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim\';
%masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-CenterStim-Fast\';
%% camera and donut corrections and FFC
%load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat','donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\ratioCorrectionFinalsTK.mat','ratioCorrectionFinal24');
ratioCorrectionFinal=ratioCorrectionFinal24;

%load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera correctionsTK 2020-03-23.mat');

%% scale bar info
umPerPix = distanceScale(60);  %=0.21 um/pixel
desiredLengthOfScale=25; % in um
numPixInScaleBar=desiredLengthOfScale/umPerPix; % convert scale length in microns to pixels
scaleHeight=12.5; %in pix
%% Create rotation matrix for converting from micron xy to pixel xy
rotationAngle=-1*-.0124;  % Multiply by negative one because we are doing the reverse transformation here
rotationMatrix=[cos(rotationAngle) -1*sin(rotationAngle); sin(rotationAngle) cos(rotationAngle)];
%% identify cellTypeStimTypeAll across all data folders
rawDataDir=getDirectories([masterRoot 'Raw Data'],'2020-07-22_PLB-Cdc42-CenterStim' )'
%rawDataDir(1)=[];

for numRawDir=1:length(rawDataDir)
    parentFolder=rawDataDir{numRawDir};
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    load([scriptRoot filesep 'keepList1.mat']);
    dataSubdir=getSubdirectories(rawDataroot,'centroid'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir);
    tempGoodCells=dataSubdir(keepList1);
    
    %define unique cell treatmetns ie KO v regular
    %cellTypeLabel=unique(getTokens(tempGoodCells,'PLB[\_\-]([\w]+\-[\w]+)'));
    cellLabelsAll=regexp(tempGoodCells,'(?<=PLB-).*(?=Cell)','match');
    cellTypeLabel=unique(cellfun(@(x) x,cellLabelsAll));
    cellTypeLabels{numRawDir}=cellTypeLabel;
    %cellTypeStimType{numRawDir}=concatenateStringsAcross2Variables(retLabelsAll,cellTypeLabel);
end

cellTypeLabelsCat=horzcat(cellTypeLabels{:});
cellTypeLabelAll=unique(cellTypeLabelsCat);
% restrict Cell type to control and C10 only


ind1=find(boolRegExp(cellTypeLabelAll,'CP-Cdc42KO'));
ind2=find(boolRegExp(cellTypeLabelAll,'SCR-SC624'));
ind3=find( boolRegExp(cellTypeLabelAll,'Pak1KO'));
badInd=[ind1 ind2 ind3];
cellTypeLabelAll(badInd)=[];



%define counter and empty image cell arrays for the block below
counter=zeros(1,length(cellTypeLabelAll));
allCount=0;


%% Specify data folder and load filenames and img Corrections

% Make FRET Before-After image pairs and analyze change in Cdc42 activity as a function of distance from the FRAP spot
% this block has been optimized for cell center stimulation.
%bins=0:3:69; %in pixles
distcomp.feature( 'LocalUseMpiexec', false );
rawDataDir=getDirectories([masterRoot 'Raw Data'],'2020-07-22_PLB-Cdc42-CenterStim' )
%rawDataDir(1)=[];
%exclude the first folder

countSubdir=0; % counter for total number of subdirectories read, allows me to build a big struct for 3 days of exp
binSize=1;
imgName1='TIRF-TomKat-1';
imgName2='TIRF-TomKat-2';

% cell number counter

 numRawDir=1;
    
    clear dataSubdir; clear keepList1;
    clear p; clear frapMask;  clear goodCells;
    parentFolder=rawDataDir{numRawDir};
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    dataSubdir=getSubdirectories(rawDataroot,'centroid'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir);
    load([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
    load([scriptRoot filesep 'frapMask1PXrt.mat'],'frapMask1PXrt','frapIndrt');
    %load([scriptRoot filesep 'movieKeepLists.mat']);
    goodCells=dataSubdir;
    
    foldDate=str2double(regexp(parentFolder, '^[0-9]+','match'));
    if foldDate < 2020
        load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat','donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');
    else
        load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera correctionsTK 2020-03-23.mat');
    end
    frapMask=frapMask1PXrt;
    frapMask=logical(frapMask);
    frapInd=frapIndrt;

   
    
    %getDataStr
    dateStr=getTokens(parentFolder,'(^[0-9]+\-[0-9]+\-[0-9]+)');
    %define counter based on cell type and frap laser power
    
    %build the metadata structure for this exp folder
    metaFrame = parseFrameMetaData(rawDataroot);
    %% define cell folder label and find it in the goodCells array
    cellName='B02_PLB-Ctrl-SC344-SC147-Ret_centroid_3mW_Cell15';
    
    kInd=find(boolRegExp(goodCells,cellName));
    
    k=kInd;
    
    % Read stage positions
    folder=[rawDataroot filesep goodCells{k}];%use dataSubdir or dataSubdirCentroid to specify all or some of your data
    % if isempty(find(boolRegExp(notKeepersFileNames, goodCells{k})))
    numImFiles=findNumImgNamesInFolder(rawDataroot, goodCells,{imgName1,imgName2});
    numBefore=findNumImgBeforeFrapStim(folder,{imgName1,imgName2},'imgbeforestim',false);
    % id the stimulus group from file name
    %pull out stim name from current folder
    stimLabel=getTokens(goodCells{k},'([0-9]+mW)');
    cellTypeLabel=regexp(goodCells{k},'(?<=PLB-).*(?=Cell)','match');
    %generate logical map for the identified stimGroup.
    cellStimLogicals=strcmpi(cellTypeLabelAll,cellTypeLabel);
    if ~isempty(find(cellStimLogicals))
        %id the position in the stimGroup logical matrix. Use this value to select
        %the correct position in res. ie res(stimInd).name;
        stimInd=find(cellStimLogicals);
        fprintf('stimInd: %i,',stimInd);
        
        
        % Read stage positions
        objMag=60;
        [pos1, pos1PassFlag]=createPos1(folder,objMag,rotationMatrix);
        
        
        % build centerMask
        % identify the cell at the center (use for frap imaging only)
        centerMask=frapMask;
        centerMask=imdilate(centerMask(:,:,1), strel('disk',50));
        
        % Read and process images
        filenames1=getFilenames(folder,imgName1);
        filenames2=getFilenames(folder,imgName2);
        cropsize=[951 921];
        %cropsize=[300 300];
        ratioImage=cell(1,length(filenames1));
        correctedRatioIm=cell(1,length(filenames1));
        grayImg=cell(1,length(filenames1));
        allOrigMasks=cell(1,length(filenames1));
        allKatIm=cell(1,length(filenames1));
        allTomIm=cell(1,length(filenames1));
        minAreaThreshTF=cell(1,length(filenames1));
        ts = getImgTimesForCellMetaData(metaFrame, filenames1,numBefore);
        if all(ts.failFlag==false) && all(diff(ts.tkTime)<4) && pos1PassFlag
            for i=1:length(filenames1)
                
                im1=double(imread([folder filesep filenames1{i}]));
                im1=(im1- medianDark1).*donutCorrection1.*halfCorrection1;
                im2=double(imread([folder filesep filenames2{i}]));
                im2=(im2-medianDark2).*donutCorrection2.*halfCorrection2;
                im = singleIms2stack(im1,im2);
                regIm=registerImagesFromQuadFit(im,p);
                cropIm2=cropImMidOut(regIm,'cropsize', cropsize);
                
                [maskFinal, cRatioIm, imGray, acceptorIm, donorIm,minAreaThreshPass]= processFretRatioImTIRF(cropIm2,centerMask,ratioCorrectionFinal,'sharpparam',[4 25],'usemultithresh',false);
                minAreaThreshTF{i}=minAreaThreshPass;
                allOrigMasks{i}=maskFinal;
                grayImg{i}=imGray;
                correctedRatioIm{i}=cRatioIm;
                allKatIm{i}=acceptorIm;
                allTomIm{i}=donorIm;
            end
            % check the minAreaThreshFilter
            minAreaFilt=cellfun(@(x) x,minAreaThreshTF);
            if all(minAreaFilt)
                
                % shift images to account for stage movement
                frapIms=makeCellArrWithDuplicateImages(frapMask,length(correctedRatioIm));
                frapIms=shiftCellImWithPos1V2(frapIms,pos1,'logical',true);
                adjFrapCoors = shiftFrapCoorsByPos1(frapInd,pos1,correctedRatioIm);
                allOrigMasks = shiftCellImWithPos1V2(allOrigMasks,pos1,'logical',true);
                grayImg=shiftCellImWithPos1V2(grayImg,pos1);
                grayImg=convertZeroPadsToNans(grayImg,allOrigMasks);
                %                     correctedRatioIm=shiftCellImWithPos1V2(correctedRatioIm,pos1);
                correctedRatioIm=convertZeroPadsToNans(correctedRatioIm,allOrigMasks);
                allKatIm=shiftCellImWithPos1V2(allKatIm,pos1);
                allKatIm=convertZeroPadsToNans(allKatIm,allOrigMasks);
                allTomIm=shiftCellImWithPos1V2(allTomIm,pos1);
                allTomIm=convertZeroPadsToNans(allTomIm,allOrigMasks);
                ratioCorrectionFinalCArray=makeCellArrWithDuplicateImages(ratioCorrectionFinal,length(correctedRatioIm));
                allRatioCorr=shiftCellImWithPos1V2(ratioCorrectionFinalCArray,pos1);
                
                
                
                % determine min and max fret range for each cell orig Range
                minFret=cellfun(@(x) prctile(vect(x),1),correctedRatioIm);
                maxFret=cellfun(@(x) prctile(vect(x),99),correctedRatioIm);
                fretRange = round([min(minFret)  max(maxFret)],2);
                
                %crop all relavant images and save. Also calc cell
                %elimination stats.
                cropsize2=[400 400];
                [imDat, failFlag] = measureCenterStimElimStatsAndSaveImDatV2(...
                    correctedRatioIm,grayImg,allOrigMasks,cropsize2,frapIms,allKatIm,allTomIm, allRatioCorr);
                if failFlag==false
                    %updat counter after script determines if the cell area is
                    %on target
                    counter(stimInd)=counter(stimInd)+1;
                    allCount=allCount+1;
                    imDat.pos1=pos1;
                    imDat.fretRange=fretRange;
                    imDat.fileName=goodCells{k};
                    imDat.parentFolder=parentFolder;
                    imDat.date=dateStr{1};
                    imDat.fileName=goodCells{k};
                    imDat.times=ts;
                    imDataSaveName=sprintf('%s-imgData.mat',goodCells{k});
                    
                    
                    % Make movie frames
                    moviePath=[saveRoot filesep goodCells{k} '.avi'];
                    xRange=[1:301];
                    yRange =[1:301];
                    tempSaveRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim\Analyzed Data\2020-07-30_PooledDatafromJun2020-july2020\Figs';
                    
                                            writeFRETMoviesFromImageSequencesPostShiftMatrix(imDat.cropFret,...
                                                imDat.cropGray, moviePath, imDat.cropfrapInd, 'imframerate',0.5,...
                                                'flip',false,'numbefore',numBefore,'objective',60,...
                                                'boundsfret',fretRange,'frapcounter',1);
                    
                    
                    foldNameFretComposite=sprintf('Fret-%s-Test',goodCells{k});
                    montagePathFretComposite=[tempSaveRoot filesep foldNameFretComposite];
                    mkdir(montagePathFretComposite)
                    
                    %define x&y range based on mask centroid
                    obJ=regionprops(imDat.cropMask{1},'Centroid');
                    cent=ceil(obJ.Centroid);
                    % function that builds a vector with a specific val in
                    % the center of function.
                    generate_array = @(CtrPt,NrPts,Res) linspace(CtrPt-Res*fix(NrPts/2), CtrPt+Res*fix(NrPts/2), NrPts);
                    NrPts = 120;
                    Res   = 1;
                   xRange= ceil(generate_array(cent(1),NrPts,Res));
                   yRange= ceil(generate_array(cent(2),NrPts,Res));
                    
                    frapSelect=[5:2:13];
                    buildCompositeFRETImgFigure(imDat.cropFret,...
                        montagePathFretComposite,frapSelect,'pos1',imDat.pos1,...
                        'frapspot',imDat.cropfrapInd,'flip',false,...
                        'numbefore',numBefore,'objective',60,...
                        'boundsfret',fretRange,'savefiletype','pdf',...
                        'imsgray',imDat.cropGray,'ts',imDat.times,'xrange',xRange,...
                        'yrange',yRange,'scalebarlength',15,'diffimgtf',true,...
                        'difframes',[1 3],'frapscattersize',50);
                   
                    
                end %if FailFlag==false
                fprintf('k=%d,', allCount);
                if mod(allCount,20)==0
                    fprintf('\n');
                end
            end %cBinOut.cellMaskPass==false
        end %if all(minAreaFilt)
    else
        fprintf('ts error\n');
    end %if all(diff(ts.tkTime)<4)
    %end %if ~isempty(find(cellStimLogicals))
    %end %if good cell is in not keeper list
    
    
    
    

%% save data

centerStimResFast20210128=uniRes;
%imgDataJunJuly112018CenterStim=uniImgDataOut;
%TempPath='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim\Analyzed Data\2020-07-30_PooledDatafromJun2020-july2020';

save([saveRoot filesep 'centerStimResFast20210128.mat'],'centerStimResFast20210128','-v7.3');


% cellsThatTiggeredUniResOut=uniResOut;
% %imgDataJunJuly112018CenterStim=uniImgDataOut;
% TempPath='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim\Analyzed Data\2020-04-24_PooledDatafromJun2018-july11-2018';
% save([TempPath filesep 'cellsThatTiggeredUniResOut.mat'],'cellsThatTiggeredUniResOut','-v7.3');
% %save('myFile.mat', 'Variablename', '-v7.3')

%% make montage from loaded imDat
%2021-04-06. I had already computed the imDat for Fret-B02_PLB-Ctrl-SC344-SC147-Ret_centroid_3mW_Cell15
% Here I drag-loaded the imDat, and I used it to make a new version of fig
% 5b.
tempSaveRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim\Analyzed Data\2020-07-30_PooledDatafromJun2020-july2020\Figs';

cellName='B02_PLB-Ctrl-SC344-SC147-Ret_centroid_3mW_Cell15';
numBefore=find(imDat.postStimInd==true,1,'first');
fretRange=[0.97 1.12];
% Make movie frames
moviePath=[tempSaveRoot filesep cellName '-test.avi'];
xRange=[1:301];
yRange =[1:301];

writeFRETMoviesFromImageSequencesPostShiftMatrix(imDat.cropFret,...
    imDat.cropGray, moviePath, imDat.cropfrapInd, 'imframerate',0.5,...
    'flip',false,'numbefore',numBefore,'objective',60,...
    'boundsfret',fretRange,'frapcounter',1);

%%
foldNameFretComposite=sprintf('Fret-%s-Test1',cellName);
montagePathFretComposite=[tempSaveRoot filesep foldNameFretComposite];
mkdir(montagePathFretComposite)

%add timesInfrapRef to imDat.times
imDat.times.timesInfrapRef=imDat.times.tkInSec-imDat.times.frapInSec;

%define x&y range based on mask centroid
obJ=regionprops(imDat.cropMask{1},'Centroid');
cent=ceil(obJ.Centroid);

% function that builds a vector with a specific val in
% the center of function.
%generate_array = @(CtrPt,NrPts,Res) linspace(CtrPt-Res*fix(NrPts/2), CtrPt+Res*fix(NrPts/2), NrPts);

NrPts = 120;
Res   = 1;
xRange= ceil(generateArray(cent(1),NrPts,Res)); % use my function instead
yRange= ceil(generateArray(cent(2),NrPts,Res));

%2021-04-06-07 I updated this buildCom... function to also calcuate fret
%ratio fold change see the varargin inputs that allow swithcing between
%differnce and fold change.
frapSelect=[5:2:13];
buildCompositeFRETImgFigure(imDat.cropFret,...
    montagePathFretComposite,frapSelect,'pos1',imDat.pos1,...
    'frapspot',imDat.cropfrapInd,'flip',false,...
    'numbefore',numBefore,'objective',60,...
    'boundsfret',fretRange,'savefiletype','pdf',...
    'imsgray',imDat.cropGray,'ts',imDat.times,'xrange',xRange,...
    'yrange',yRange,'scalebarlength',15,'diffimgtf',false,...
    'relchangetf',true,'difframes',[1 3],'frapscattersize',50,...
    'scalebartexttf', false,'rcrange',[0.91 1.09]);
%%
%% save a copy of this script in the scriptRoot
scriptRoot=montagePathFretComposite;
Version='20210407Fig5_v2';
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