%% Script for computing front-back profile of FRET sensor data for polarized cells


%2021-04-04 Update: For publication I need to submit the source data that
%was used to generate each plot. The data will be stored in an excel sheet.
%at the bottom of this script I added the section for this data export.


%Note: this code was used to generate the data and plots in Figure 1c of
%Bell et al 2021.
%%
%root = 'E:\Users\collinsLab\Documents\GB_Data\2017-11-15-PLB-PP1TK-UnderAg';  % This should be the path to the folder containing all the data from the experiment
%% set roots
masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\SC147-FrontBack\';
%% camera and donut corrections
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat',...
    'donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\ratioCorrectionFinalsTK.mat','ratioCorrectionFinal24');
%%
rawDataDir=getDirectories([masterRoot 'Raw Data'],'11-14')
countSubdir=0; % counter for total number of subdirectories read, allows me to build a big struct for 3 days of exp
binSize=1;
imgName1='TIRF-TomKat-1';
imgName2='TIRF-TomKat-2';
numRawDir=1;
rawDataDir{numRawDir};

clear dataSubdir;
parentFolder=rawDataDir{numRawDir};
% mkdir([masterRoot 'Scripts' filesep parentFolder]); % makes a folder for the exp day in the scripts folder
% mkdir([masterRoot 'Analyzed Data' filesep parentFolder]);% makes a folder for the exp day in the Analyzed Data folder
rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
dataSubdir=getSubdirectories(rawDataroot,'SC147'); % get the subdirectories for each folder in rawData
dataSubdir=natsortfiles(dataSubdir); % sort the dat
% load p,ratioCorrection
load([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
%load([scriptRoot filesep '2D ratio correction-60xTK.mat'], 'ratioCorrectionFinal');

%% Test ratio image
subNum=2;
frameNum=10;
clear im
folder=[rawDataroot filesep dataSubdir{subNum}];
filenames1=getFilenames(folder,imgName1);
filenames2=getFilenames(folder,imgName2);
im(:,:,1)=double(imread([folder filesep filenames1{frameNum}]));
im(:,:,1)=(im(:,:,1)-medianDark1).*donutCorrection1.*halfCorrection1;
im(:,:,2)=double(imread([folder filesep filenames2{frameNum}]));
im(:,:,2)=(im(:,:,2)-medianDark2).*donutCorrection1.*halfCorrection2;
im=registerImagesFromQuadFit(im,p);
im=im(50:970,30:980,:);
cropsize=[951 921];
cropIm2 = cropImMidOut(im,'cropsize', cropsize);
% background subtraction

imSum=sum(cropIm2,3);
sharpParam=[4 50];

[maskBGinitial,objectRectangles]=fretGetInitialObjectsAndBGmask_Noisy_BG(imSum);
imForSeg0=imageSubtractBackgroundWithObjectMaskingAndPadding(imSum,100,50,maskBGinitial);
[maskFinal,coors]=fretGetCellMasks_63xWithErode(imForSeg0,objectRectangles,2000, imSum,imSum,sharpParam);
CFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(cropIm2(:,:,2),100,50,maskBGinitial);
YFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(cropIm2(:,:,1),100,50,maskBGinitial);

%smoothing
YFPsmooth=YFPFinal;
CFPsmooth=CFPFinal;
YFPsmooth(~maskFinal)=nan;% mask background from cell
CFPsmooth(~maskFinal)=nan;
YFPsmooth=ndnanfilter(YFPsmooth,fspecial('gaussian',7,2),7);
CFPsmooth=ndnanfilter(CFPsmooth,fspecial('gaussian',7,2),7);
YFPsmooth(~maskFinal)=nan; %mask cell from background
CFPsmooth(~maskFinal)=nan;
% ratio image
ratioImage=YFPsmooth./CFPsmooth;
ratioImage=ratioImage./ratioCorrectionFinal;
%ratioImage(~isnan(ratioImage) & ratioImage<1.51)=1.51;
imagesc(ratioImage);
%imagesc(ratioImage,[1.5 4.5]);
colormap(parula_black);
[nansum(YFPsmooth(:))/nansum(CFPsmooth(:)) nanmean(ratioImage(:)) nanmedian(ratioImage(:))]
%% Generate Mask and Coordinates for each image

rawDataDir=getDirectories([masterRoot 'Raw Data'],'11-14')
countSubdir=0; % counter for total number of subdirectories read, allows me to build a big struct for 3 days of exp
binSize=1;
imgName1='TIRF-TomKat-1';
imgName2='TIRF-TomKat-2';
for numRawDir=1:length(rawDataDir)
    clear dataSubdir; clear p; clear ratioCorrectionFinal; clear fretData
    parentFolder=rawDataDir{numRawDir};
    mkdir([masterRoot 'Scripts' filesep parentFolder]); % makes a folder for the exp day in the scripts folder
    mkdir([masterRoot 'Analyzed Data' filesep parentFolder]);% makes a folder for the exp day in the Analyzed Data folder
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    dataSubdir=getSubdirectories(rawDataroot,'SC147');  % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the dat
    % load p,ratioCorrection
    load([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
    load([scriptRoot filesep 'keepList1.mat'],'keepList1');
    clear allCoors; clear allTraj; clear keepList;
    cropSize=[600 600];
    ratioCorrectionFinal=cropImMidOut(ratioCorrectionFinal24,'cropsize',cropSize);
    
    
   % remove dataSubdir with no images 
   
   if exist('keepList1')==1
       keepList=keepList1;
   else
       
       for t=1:length(dataSubdir)
           folder=[rawDataroot filesep dataSubdir{t}];
           if any(size(dir([folder '/*.tif' ]),1)) ==1
               filenames1=getFilenames(folder,imgName1);
               filenames2=getFilenames(folder,imgName2);
               if ~isempty(filenames1) & length(filenames1) >10
                   keepList2(t)=1;
               end
           end
           
       end
       keepList=find(keepList2);
   end
    
   dataSubdir=dataSubdir(keepList);
   
    for i=1:length(dataSubdir)
        folder=[rawDataroot filesep dataSubdir{i}];
        clear maskBGinitial;
        clear maskFinal;
        clear coors;
        clear grayImages; clear fretImages; ;
        bgBlockSize=100;
        bgPadSize=50;
        filenames1=getFilenames(folder,imgName1);
        filenames2=getFilenames(folder,imgName2);
       
        cropSize=[600 600]; sharpParam=[4 50];
        fretImages=cell(1,length(filenames1));
        grayImages=cell(1,length(filenames1));
        coors=cell(1,length(filenames1));
        maskFinal=cell(1,length(filenames1));
        parfor t=1:length(filenames1)
           
            im1=double(imread([folder filesep filenames1{t}]));
            im1=(im1- medianDark1).*donutCorrection1.*halfCorrection1;
            im2=double(imread([folder filesep filenames2{t}]));
            im2=(im2-medianDark2).*donutCorrection2.*halfCorrection2;
            im = singleIms2stack(im1,im2);
            regIm=registerImagesFromQuadFit(im,p);
            cropIm=regIm(50:970,30:980,:);
            %Crop again
            
            cropIm2 = cropImMidOut(cropIm,'cropsize', cropSize);
            
            
            % background subtraction
            imSum=sum(cropIm2,3);
            try
                [maskBGinitial,objectRectangles]=fretGetInitialObjectsAndBGmask(imSum);
            catch
                objectRectangles={[1 1 601 601]};
                
            end
            if ~isempty(objectRectangles)
                
                imForSeg0=imageSubtractBackgroundWithObjectMaskingAndPadding(imSum,150,100,maskBGinitial);
                %[maskFinal,~]=fretGetCellMasks(imForSeg0,{[1 1 1024 1024]},100);
                [ maskF,coors{t}]=fretGetCellMasks_63xWithErode(imForSeg0,objectRectangles,2500, imSum,imSum,sharpParam);
                sizeMask=regionprops( maskF,'Area'); % remove images with no cells. In these cases, the channel becomes the whole mask.
                maskInd=arrayfun(@(x) (x.Area>20000), sizeMask);
                
                if sum(maskInd)<1
                    
                    maskFinal{t}=bwareaopen( maskF,2000);
                    
                    CFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(cropIm2(:,:,2),150,100,maskBGinitial);
                    YFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(cropIm2(:,:,1),150,100,maskBGinitial);
                    
                    %smoothing
                    YFPsmooth=YFPFinal;
                    CFPsmooth=CFPFinal;
                    YFPsmooth(~maskFinal{t})=nan;% mask background from cell
                    CFPsmooth(~maskFinal{t})=nan;
                    YFPsmooth=ndnanfilter(YFPsmooth,fspecial('gaussian',7,2),7);
                    CFPsmooth=ndnanfilter(CFPsmooth,fspecial('gaussian',7,2),7);
                    YFPsmooth(~maskFinal{t})=nan; %mask cell from background
                    CFPsmooth(~maskFinal{t})=nan;
                    % ratio image
                    ratioImage=YFPsmooth./CFPsmooth;
                    fretImages{t}=ratioImage./ratioCorrectionFinal;
                    grayImages{t}=imForSeg0;
                    %not used anymore           fretData(i).meanWeighted(t)=nansum(YFPsmooth(:))/nansum(CFPsmooth(:));
                    
                else
                    fretImages{t}=nan(size(cropIm2,1), size(cropIm2,2));
                grayImages{t}=nan(size(cropIm2,1), size(cropIm2,2));
                maskFinal{t}=nan(size(cropIm2,1), size(cropIm2,2));
                end %if sumMask1
                
            else
                fretImages{t}=nan(size(cropIm2,1), size(cropIm2,2));
                grayImages{t}=nan(size(cropIm2,1), size(cropIm2,2));
                maskFinal{t}=nan(size(cropIm2,1), size(cropIm2,2));
            end %if empty objectRect
           
                
        end %parfor
       
        maxVals=max(cellfun(@(x) prctile(vect(x),99),fretImages));
        minVals=min(cellfun(@(x) prctile(vect(x),1),fretImages));
        FretRange=[mean(minVals), mean(maxVals)];
        %writing movie
        %fbDataAll{i}{r}.yOut(frameCount,:)=yOut;
        filePath=[saveRoot filesep dataSubdir{i} '_movie.avi'];
        writeFRETMoviesFromImageSequences(fretImages,grayImages,filePath,'boundsfret',FretRange,'objective',60,...
            'imframerate',3,'numbefore',5);
        %save([saveRoot filesep 'masks and coors.mat'], 'maskFinal', 'coors', 'maskBGinitial');
        allCoors{i}=coors; %collating coordinates
        allTraj{i}=basicTrackingFunction(allCoors{i},100); %tracking cells
        fretData(i).mean=cellfun(@(x) nanmean(x(:)),fretImages);
        fretData(i).median=cellfun(@(x) nanmedian(x(:)),fretImages);
        %fretData(i).grayImages=grayImages;
        fretData(i).fretImages=fretImages;
        fretData(i).maskFinal=maskFinal;
        fretData(i).fretRange=FretRange;
        fretData(i).dataSubdir=dataSubdir{i};
 
        
        fprintf('f=%i,',i);
    end %for 1;length(dataSubdir)
    save([saveRoot filesep 'fretData Weighted Mean.mat'], 'fretData','-v7.3');
    save([saveRoot filesep 'allCoors and allTraj.mat'], 'allCoors', 'allTraj');
end %numRawDir

%% Load fret data and allCoors/allTraj
 clear fretData allCoors allTraj dataSubdir
 rawDataDir=getDirectories([masterRoot 'Raw Data'],'FrontBack')
 
for numRawDir=1:length(rawDataDir)
    parentFolder=rawDataDir{numRawDir};
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    load([saveRoot filesep 'fretData Weighted Mean.mat'], 'fretData');
    load([saveRoot filesep 'allCoors and allTraj.mat'], 'allCoors', 'allTraj');
    fretDat(numRawDir).fretData=fretData;
    coorNtraj(numRawDir).allCoors=allCoors;
    coorNtraj(numRawDir).allTraj=allTraj;
    dataSubdir=getSubdirectories(rawDataroot,'SC147'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the dat
    allDataSubDir(numRawDir).dataSubdir=dataSubdir;
    %clear fretData allCoors allTraj
end

fretData=[fretDat(1).fretData fretDat(2).fretData fretDat(3).fretData];
dataSubdir=[allDataSubDir(1).dataSubdir allDataSubDir(2).dataSubdir allDataSubDir(3).dataSubdir];
allCoors=[coorNtraj(1).allCoors coorNtraj(2).allCoors coorNtraj(3).allCoors];
allTraj=[coorNtraj(1).allTraj coorNtraj(2).allTraj coorNtraj(3).allTraj];

%% Save combo fretData amd AllCOors
 save([saveRoot filesep 'Combo FretData.mat'], 'fretData','-v7.3');
    save([saveRoot filesep 'Combo allCoors and allTraj.mat'], 'allCoors', 'allTraj');
%% Load data
 load([saveRoot filesep 'Combo FretData.mat'], 'fretData');
    load([saveRoot filesep 'Combo allCoors and allTraj.mat'], 'allCoors', 'allTraj');

%% Define experimental groups
groups={[1:length(fretData)]};

%% Pick cells
minTrajLen=5;

minFrame=2;
maxFrame=59;
recordsByWell=cell(size(groups));
wells=sort(vect(cell2mat(groups))');

for i=1:length(wells)
    recordsByWell{i}=[];
    numFrames=length(allCoors{i});
    for j=1:length(allTraj{i})
        if sum(allTraj{i}{j}(:,end)>=minFrame)>=minTrajLen
            % define frame first and then check for min area
            frameFirst=max(minFrame,allTraj{i}{j}(1,end));
            
            if allTraj{i}{j}(1,3) >4000 
            figure(1);
            set(gcf,'Position',[50 280 1300 700]);
            
            % Cell area vs time
            subplot('Position',[0.05 0.77 0.35 0.18]);
            plot(allTraj{i}{j}(:,end),allTraj{i}{j}(:,3)); title('Cell Area');
            xlim([minFrame-2, numFrames]); ylim([0 2e4]); set(gca,'box','off');
            
            % Cell speed vs time
            speedV=sqrt( (allTraj{i}{j}(3:end,1) - allTraj{i}{j}(1:(end-2),1)).^2 + (allTraj{i}{j}(3:end,2) - allTraj{i}{j}(1:(end-2),2)).^2 );
            speedV=smooth(speedV,7,'sgolay');
            subplot('Position',[0.05 0.53 0.35 0.18]);
            plot(allTraj{i}{j}(2:(end-1),end),speedV,'r'); title('Cell Speed');
            xlim([minFrame-2 numFrames]); set(gca,'box','off');
            
            
            % Cell position
            subplot('Position',[0.05 0.29 0.35 0.18]);
            plot(allTraj{i}{j}(:,end),allTraj{i}{j}(:,1)); title('X Position');
            xlim([minFrame-2 numFrames]); set(gca,'box','off');
            subplot('Position',[0.05 0.05 0.35 0.18]);
            plot(allTraj{i}{j}(:,end),allTraj{i}{j}(:,2)); title('Y Position');
            xlim([minFrame-2 numFrames]); set(gca,'box','off');
            set(gca,'YDir','reverse');
            
            % Images for reference
     
            positions={[0.415 0.55 0.13 0.4],[0.565 0.55 0.13 0.4],[0.715 0.55 0.13 0.4],[0.865 0.55 0.13 0.4], ...
                [0.415 0.05 0.13 0.4],[0.565 0.05 0.13 0.4],[0.715 0.05 0.13 0.4],[0.865 0.05 0.13 0.4]};
            subplot('Position',positions{1});
             frameFirst=max(minFrame,allTraj{i}{j}(1,end));
            imagesc(fretData(i).fretImages{frameFirst});
            set(gca,'XTick',[]); set(gca,'YTick',[]);
            if ~isempty(allTraj{i}{j}(find(allTraj{i}{j}(:,end)==frameFirst,1),4:7));
                rectangle('Position',allTraj{i}{j}(find(allTraj{i}{j}(:,end)==frameFirst,1),4:7),'EdgeColor',[1 1 0]);
                title(['Frame ' num2str(frameFirst)]);
                
                subplot('Position',positions{end});
                frameLast=min(maxFrame,allTraj{i}{j}(end,end));
                imagesc(fretData(i).fretImages{frameLast});
                set(gca,'XTick',[]); set(gca,'YTick',[]);
                rectangle('Position',allTraj{i}{j}(find(allTraj{i}{j}(:,end)==frameLast,1),4:7),'EdgeColor',[1 1 0]);
                title(['Frame ' num2str(frameLast)]);
                
                for k=1:6
                    subplot('Position',positions{k+1});
                    frameMid=round(((7-k)*frameFirst+k*frameLast)/7);
                    imagesc(fretData(i).fretImages{frameMid});
                    set(gca,'XTick',[]); set(gca,'YTick',[]);
                    rectangle('Position',allTraj{i}{j}(find(allTraj{i}{j}(:,end)==frameMid,1),4:7),'EdgeColor',[1 1 0]);
                    title(['Frame ' num2str(frameMid)]);
                end
            end %if is empty rectangle traj
            fprintf('Well = %i  Cell = %i\n',i,j);
            temp='x';
            while (~isempty(temp) && ~boolRegExp({temp},'[0-9]+ [0-9]+'))
                %Frames: startframe endframe convention.
                temp=input('Frames: ','s');
            end
            if ~isempty(temp)
                temp=str2double(getTokens(temp,'([0-9]+)'));
                recordsByWell{i}=[recordsByWell{i}; i j temp];
            end
            close(1);
            end %if area >4000
        end %if minFRame is > frame length
    end %j
end % i
%% Organize cell choice records % this is working 2/19/16 gb
records=cell(size(groups));
for i=1:length(groups)
    records{i}=cell2mat(vect(recordsByWell(groups{i})));
end
%% keep the older records for safe keeping
oldRecords=records;
%% screen records by well to determine if there are frames with no mask
folderNum=0; imageNum=0;
records=oldRecords;

for i=1:length(records)
    for r=1:size(records{i},1)
        if folderNum~=records{i}(r,1)
            folderNum=records{i}(r,1);
            imageNum=0; %So that I know it's a new image
        end
       clear maskCheck fretCheck;
       if records{i}(r,4)<records{i}(r,3)
           records{i}(r,4)=records{i}(r,3)+1;
       end
       
       
       for imageNum=records{i}(r,3):records{i}(r,4)
           imNumCount=imNumCount+1;
           fretCheck(imageNum)=~isempty(find(any(~isnan(fretData(folderNum).fretImages{imageNum}))));
           maskCheck(imageNum)=~isempty(find(any(~isnan(fretData(folderNum).maskFinal{imageNum}))));
       end
           logCheck=isempty(find(fretCheck~=maskCheck));
           if ~logCheck
               fretConsec=max(numConsecutiveInList(fretCheck));
               maskConsec=max(numConsecutiveInList(maskCheck));
               if maskConsec<fretConsec
                   masterCheck=maskCheck;
               else 
                   masterCheck=fretCheck;
               end
           else
               masterCheck=maskCheck;
           end
           out=numConsecutiveInList(masterCheck);
           sizeConsec=out(out>0);
           largestConsec=max(sizeConsec);
           ind=find(out==largestConsec);
           ind=ind(1);
           imgKeep=[ind:(ind+largestConsec)-2]; %-1 for counting -1 again because distmask requires the frame after to be nan free too.
           records{i}(r,3)=imgKeep(1);
           records{i}(r,4)=imgKeep(end);
           
           
        
    end
end
recordsTopCopy=records;

%% Save cell choice records
save([saveRoot filesep 'FrontBackAnalysisRecordsTopCopy.mat'],'records','recordsTopCopy','recordsByWell');
%% Load cell choice records
load([saveRoot filesep 'FrontBackAnalysisRecords.mat'],'records','recordsByWell');
%%
load([saveRoot filesep 'FrontBackAnalysisRecordsTopCopy.mat'],'records','recordsTopCopy','recordsByWell');

%% Load combo fretData amd AllCOors
 load([saveRoot filesep 'Combo FretData.mat'], 'fretData');
    load([saveRoot filesep 'Combo allCoors and allTraj.mat'], 'allCoors', 'allTraj');
%% Compute CDC42 intracellular gradient data
smoothFilt=fspecial('disk',2.5)>0.015; smoothFilt=smoothFilt/sum(smoothFilt(:));

%the lenght of bins seems too long. will restrict this to ~69 bins
bins=0:2.5:300;
%bins=0:3:69;
fixedN=200;

folderNum=0; imageNum=0;
figure(1); set(1,'Position',[50 200 400 600]);
figure(2); set(2,'Position',[500 200 500 500]);
for i=1:length(records)
    for r=1:size(records{i},1)
        if folderNum~=records{i}(r,1)
            folderNum=records{i}(r,1);
            imageNum=0; %So that I know it's a new image
%             folder=[rawDataroot filesep dataSubdir{folderNum}];
%             filenames1=getFilenames(folder,'TIRF-TomKat-1');
%             filenames2=getFilenames(folder,'TIRF-TomKat-2');
%             load([saveRoot filesep 'masks and coors.mat'],'maskBGinitial','maskFinal','coors');
        end
        frameCount=0;
        for imageNum=records{i}(r,3):records{i}(r,4)
            clear maskFinal;
            %frameCount=frameCount+1;
            if imageNum < length(fretData(folderNum).maskFinal) & imageNum >1
                figure(1); imagesc(fretData(folderNum).fretImages{imageNum}); 
                %TITLE the displayed image
                titleStr=sprintf('%s imgNum %s',fretData(folderNum).dataSubdir, num2str(imageNum));
                title(regexprep(titleStr,'_',' '));
                %remove small objects that jack up the mask
                maskFinal=bwareaopen(fretData(folderNum).maskFinal{imageNum},300);
                objects=regionprops(maskFinal,'PixelIdxList','PixelList','Centroid','BoundingBox');
                thisTraj=allTraj{records{i}(r,1)}{records{i}(r,2)};
                if ~isempty(find(thisTraj(:,end)==imageNum,1))
                    cellCent=round(thisTraj(find(thisTraj(:,end)==imageNum,1),1:2));
                    cellNum=find(round(arrayfun(@(x) x.Centroid(1),objects))==cellCent(1) & round(arrayfun(@(x) x.Centroid(2),objects))==cellCent(2));
                    if ~isempty(objects(cellNum))
                        frameCount=frameCount+1;
                    cx1=floor(max(1,objects(cellNum).BoundingBox(1)));
                    cy1=floor(max(1,objects(cellNum).BoundingBox(2)));
                    cx2=ceil(min(size(maskFinal,2),sum(objects(cellNum).BoundingBox([1 3]))));
                    cy2=ceil(min(size(maskFinal,1),sum(objects(cellNum).BoundingBox([2 4]))));
                    thisMask=false(size(maskFinal));
                    thisMask(objects(cellNum).PixelIdxList)=true;
                    pixInd=objects(cellNum).PixelIdxList;
                    distMask=plbMasks2DistMask(fretData(folderNum).maskFinal,imageNum); distMask(~thisMask)=nan;
                    [yOut,nOut] = evalBySlidingBins(vect(distMask(cy1:cy2,cx1:cx2)),vect(fretData(folderNum).fretImages{imageNum}(cy1:cy2,cx1:cx2)),bins,1.5*(bins(2)-bins(1)),'nanmean');
                   [xOutN,yOutN] = evalBySlidingBinsFixedN(vect(distMask(cy1:cy2,cx1:cx2)),vect(fretData(folderNum).fretImages{imageNum}(cy1:cy2,cx1:cx2)),200,'nanmean',200);
                    figure(2); hold on; plot(xOutN,yOutN,'.-','Color',[1 0 0]);
                    fbDataAll{i}{r}.yOut(frameCount,:)=yOut;
                    fbDataAll{i}{r}.nOut(frameCount,:)=nOut;
                    fbDataAll{i}{r}.xOutN{frameCount}=xOutN;
                    fbDataAll{i}{r}.yOutN{frameCount}=yOutN;
                    end %if ~ispemtpy objects(cell num)
                end % if is empty traj for cell center
            end% if imageNum exceeds length of images
        end
    end
end
%% save FBData all
save([saveRoot filesep 'FBDataAll.mat'],'fbDataAll');

%% load FBData
load([saveRoot filesep 'FBDataAll.mat'],'fbDataAll');

%% filter yOut values use for bins=0:2.5:300;
%The large bin size means that we first need to deal with nan filled bins


fbDataAllGood=fbDataAll;
i=1; % for cell line
for r=1:size(records{i},1) % for each cell
    sizeyOut=size(fbDataAllGood{1}{r}.yOut,1); % number of frames
    for s=1:sizeyOut  % for each frame
      lastValPos=find(~isnan(fbDataAllGood{1}{r}.yOut(s,:)),1,'last');
      if isempty(lastValPos)
           keepers(1,r,s)=0;
      else
          %check to makesure its not all zeros
          if ~isempty(find(fbDataAllGood{1}{r}.yOut(s,:)>0))
           %lastValPos=38;
           % only checks last val
%            lastVal=fbDataAllGood{1}{r}.yOut(s,lastValPos);
%            firstVal=fbDataAllGood{1}{r}.yOut(s,1);
%            keepers(1,r,s)=firstVal>lastVal;
      % try for any vals greater than val1
       lastVal=fbDataAllGood{1}{r}.yOut(s,lastValPos);
       lastValAll(i,r,s)=lastVal;
       notFirst=fbDataAllGood{1}{r}.yOut(s,2:lastValPos); % every value except first
           firstVal=fbDataAllGood{1}{r}.yOut(s,1); % first value
           % implement negative spike filter for cell rear
           
           % if the FRET-ratio drops below the expected values 
           if sum(notFirst<0.85)>0 
           tempDiff = diff(notFirst(35:end)); % ++++++++++++++++++ TO-DO: check this value for accuracy later
           dipIndex = 34+find(tempDiff == min(tempDiff(:))); % find the value within the window where steepest drop occurs
           end % if sum(notFirst<...
           
           sumKeep=sum(firstVal<notFirst); % discard any cells where first value is not max
           keepers(1,r,s)=sumKeep<1; 
      
          end %if check for fretVals
      end %if not all nans
      if keepers(1,r,s)
          fbDataFinal{1}{r}.yOut(s,:)=fbDataAllGood{1}{r}.yOut(s,:);
          fbDataFinal{1}{r}.nOut(s,:)=fbDataAllGood{1}{r}.nOut(s,:);
          recInd(r)=r;
          
      end
    end
end


%% filter yOut values where the fret ratio for yOut(1) < yOut(end)
%use this block for when bins are set to bins=0:3:69;

fbDataAllGood=fbDataAll;
i=1;
for r=1:size(records{i},1)
    yGood=find(fbDataAllGood{1}{r}.yOut(:,1)>fbDataAllGood{1}{r}.yOut(:,end));
    fbDataAllGood{1}{r}.yOut=fbDataAllGood{1}{r}.yOut(yGood,:);
    fbDataAllGood{1}{r}.nOut=fbDataAllGood{1}{r}.nOut(yGood,:);
%     yGood2=find(fbDataAllGood{1}{r}.yOut(:,1)>1.1);
%     fbDataAllGood{1}{r}.yOut=fbDataAllGood{1}{r}.yOut(yGood2,:);
%     fbDataAllGood{1}{r}.nOut=fbDataAllGood{1}{r}.nOut(yGood2,:);
    
    sizeyOutN=length(fbDataAllGood{1}{r}.yOutN);
    for s=1:sizeyOutN
        clear yGoodN
        yGoodN=fbDataAllGood{1}{r}.yOutN{s}(1)>fbDataAllGood{1}{r}.yOutN{s}(end);
        if yGoodN
            fbDataAllGood{1}{r}.yOutN{s}=fbDataAllGood{1}{r}.yOutN{s};
            fbDataAllGood{1}{r}.xOutN{s}=fbDataAllGood{1}{r}.xOutN{s};
        else
            fbDataAllGood{1}{r}.yOutN{s}=[];
            fbDataAllGood{1}{r}.xOutN{s}=[];
        end
    end
end


%%
save([saveRoot 'front-back gradients_' datestr(now,'yyyy-mm-dd') '.mat'],'fbDataAll','fbDataAllGood','fbDataFinal','bins','fixedN','-v7.3');
%%
load([root 'front-back gradients_2013-05-30.mat'],'fbDataAll','bins','fixedN');
%% Add two more fields to the data
interpBins=1:2:150;
ptiles=2:2:98;
for i=1:length(fbDataAll)
    for j=1:length(fbDataAll{i})
        interpMat=[];
        ptileMat=[];
        for k=1:length(fbDataAll{i}{j}.xOutN)
            interpMat(k,:)=interp1(fbDataAll{i}{j}.xOutN{k},fbDataAll{i}{j}.yOutN{k},interpBins);
            ptileMat(k,:)=interp1(100*(1:length(fbDataAll{i}{j}.xOutN{k}))/length(fbDataAll{i}{j}.xOutN{k}),fbDataAll{i}{j}.yOutN{k},ptiles);
        end
        fbDataAll{i}{j}.interpBins=interpBins;
        fbDataAll{i}{j}.interpMat=interpMat;
        fbDataAll{i}{j}.ptiles=ptiles;
        fbDataAll{i}{j}.ptileMat=ptileMat;
    end
end

%% Distance Scale
umPerPix = distanceScale(60);  %=0.21 um/pixel for 60x
     
%% Plot all traces for one cell
figure; hold on;
groupNum=1;
cellNum=46;
for k=1:size(fbDataFinal{groupNum}{cellNum}.nOut,1)
    plot(bins*umPerPix,fbDataFinal{groupNum}{cellNum}.yOut(k,:));
end
%%
%% Plot all cell-averaged traces for one or more conditions
%clear fbDataFinal;
%fbDataFinal=fbDataAll;
figure; hold on;
labels={'SC147'};
colors=parula(54);
bins=0:2.5:300;
x=bins*umPerPix;
i=1;
%lineCounter=
for r=1:size(records{i},1)
    if ~isempty(fbDataFinal{i}{r})
    v=nanmean(fbDataFinal{i}{r}.yOut,1);
   % yyaxis left
    plot(x(1:34),v(1:34)./mean(v(1:5)));
%     dx = mean(diff(x));                                 % Find Mean Differece In ‘x’ Values
%     dy = gradient(v./mean(v(1:5)),dx);
%     %yyaxis right
    %plot(x,dy)
    %         v=nanmean(fbDataAll{groupNum}{j}.ptileMat,1);
    %         plot(ptiles,v,'Color',colors{groupNum});
    end
end

xlim([0 18])

dx = mean(diff(x));                                 % Find Mean Differece In ‘x’ Values
dy = gradient(y,dx);                                % Calculate Slope Of Data
xq = find((x >= 0.275) & (x <= 0.325));             % Index Of Area Of Interest ‘x’ Values
slope_x = x(xq);                                    % ‘Slope’ ‘x’ Values Of Interest
slope_y = dy(xq);                                   % ‘Slope’ ‘y’ Values Of Interest
figure(2)
plot(x, dy, '-g')
hold on
plot(slope_x, slope_y, '-r', 'LineWidth',1)
hold off
grid
legend('Slope Of All Data', 'Slope Of Data In Region Of Interest', 'Location','W')

%% Plot all cell-averaged traces for one or more conditions
%clear fbDataFinal;
%fbDataFinal=fbDataAll;
figure; hold on;

colors=parula(length(fbDataFinal{1}));
bins=0:2.5:300;
i=1;
%lineCounter=
xlim([0 30])
ylim([0.85 1.05])
for r=1:10%size(records{i},1)

    
    if ~isempty(fbDataFinal{i}{r})
    v=nanmean(fbDataFinal{i}{r}.yOut,1);
   
    plot(bins*umPerPix,v./mean(v(1:5)));
    labels{r}=sprintf('%s',num2str(r));
    %         v=nanmean(fbDataAll{groupNum}{j}.ptileMat,1);
    %         plot(ptiles,v,'Color',colors{groupNum});
    %pause
    end
end

legend(labels)
title('Cdc42 Activity v Distance from Protrusion', 'FontSize', 16)
ylabel('Mean FRET Ratio');
xlabel('Distance from Protrusion (\muM)')
%%  Plot averaged results using filtered variable fbDataFinal
bins=0:2.5:300;
binsInuM=bins*umPerPix; % truncate at 20 um.
g=figure; hold on;
labels={'SC147'};
colors={[0 0.5 0.75]};
groupNum=1;
    vmat=nan(length(fbDataFinal{groupNum}),length(binsInuM));
    temp=nan(length(fbDataFinal{groupNum}),length(binsInuM));
    for j=1:length(fbDataFinal{groupNum})
        clear temp
        if ~isempty(fbDataFinal{groupNum}{j})
            temp=nanmean(fbDataFinal{groupNum}{j}.yOut,1);
        vmat(j,:)=temp/nanmean(temp(1:5));
        %temp(res(stimG).numExp<thresh)=nan; %change 5 to 2
        end
    end
    %nans kill the patch function inside drawShadedErrorRegion. Need to
    %remove nans from nanmeanVmat
    meanVmat=nanmean(vmat,1);
    vmatInd=find(~isnan(meanVmat),1,'last');
    meanVmat=meanVmat(1:vmatInd);
    errvals=nanstd(vmat,0,1)./sqrt(sum(~isnan(vmat),1));
    errvals=errvals(1:vmatInd);
    drawShadedErrorRegion(binsInuM(1:vmatInd),meanVmat,errvals,colors{groupNum});
   h=plot(binsInuM(1:vmatInd),meanVmat,'-','Color',colors{groupNum}, 'LineWidth',2);

legend(labels);
xlim([0 20]);
%ylim([0.97 1.03])
title('Mean Cdc42 Activity v Distance from Protrusion', 'FontSize', 16)
ylabel('Mean FRET Ratio');
xlabel('Distance from Protrusion (\muM)')
%% Plot averaged results
g=figure; hold on;
labels={'SC62'};
colors={[0 0.5 0.75]};
groupNum=1;
    vmat=nan(length(fbDataAllGood{groupNum}),length(bins));
    temp=nan(length(fbDataAllGood{groupNum}),length(bins));
    for j=1:length(fbDataAllGood{groupNum})
        vmat(j,:)=nanmean(fbDataAllGood{groupNum}{j}.yOut,1);
        %temp(res(stimG).numExp<thresh)=nan; %change 5 to 2
        
    end
    errvals=nanstd(vmat,0,1)./sqrt(sum(~isnan(vmat),1));
    drawShadedErrorRegion(bins*umPerPix,nanmean(vmat,1),errvals,colors{groupNum});
   h=plot(bins*umPerPix,nanmean(vmat,1),'-','Color',colors{groupNum}, 'LineWidth',2);

legend(labels);
xlim([0 20]);
ylim([0.9 1.35])
title('Mean Cdc42 Activity v Distance from Protrusion', 'FontSize', 16)
ylabel('Mean FRET Ratio');
xlabel('Distance from Protrusion (\muM)')
%% %%
vmat147=vmat;
errvals147=errvals;
binsInuM147=binsInuM;
%% save figure and graph data
save([saveRoot filesep 'SC147 FrontBack graph data'],'vmat147','errvals147','binsInuM147');
fileNameMat=sprintf('SC147 FrontBack Plot.fig');
fileNamePDF=sprintf('SC147 FrontBack Plot.PDF');
savefig(g,[saveRoot filesep fileNameMat]);
printFigurePDF([saveRoot filesep fileNamePDF],'-bestfit');


%% load SC62 data
loadRoot='E:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\Cdc42-SC62-FrontBack\Analyzed Data\2019-10-31-PLB-SC62-FrontBack';
load([loadRoot filesep 'SC62 FrontBack graph data'],'vmat62','errvals62','binsInuM62');

%%
vmat62=vmat;
errvals62=errvals;
binsInuM62=binsInuM;
%% Plot averaged results
g=figure; hold on;
labels={'SC62'};
colors={[0 0.5 0.75]};
groupNum=1;
    
    drawShadedErrorRegion(binsInuM62,nanmean(vmat62,1),errvals62,colors{groupNum});
   h=plot(binsInuM62,nanmean(vmat62,1),'-','Color',colors{groupNum}, 'LineWidth',2);

legend(labels);
xlim([0 20]);
ylim([0.65 1.1])
title('Mean Cdc42 Activity v Distance from Protrusion', 'FontSize', 16)
ylabel('Mean FRET Ratio');
xlabel('Distance from Protrusion (\muM)')

%% save figure and graph data
save([saveRoot filesep 'SC62 FrontBack graph data'],'vmat62','errvals62','binsInuM62');
fileNameMat=sprintf('SC62 FrontBack Plot.fig');
fileNamePDF=sprintf('SC62 FrontBack Plot.PDF');
savefig(g,[saveRoot filesep fileNameMat]);
printFigurePDF([saveRoot filesep fileNamePDF],'-bestfit');

%% plot SC62 adn Sc147
colors={[0 0.5 0.75; 0    0.5000    0.4000]};
colororder=colors;
%errvalsSC147=errvals;
sc147Y=nanmean(vmat147,1);
sc147Y=sc147Y./sc147Y(1); % GB added 8/30/21 to check normalization
g=figure;
yyaxis left
drawShadedErrorRegion(binsInuM147,sc147Y,errvals147,colors{1}(1,:));
plot(binsInuM147,sc147Y,'-','Color',colors{1}(1,:), 'LineWidth',2);
ylabel('TomKat Sensor Mean FRET Ratio');
ylim([.9 1]);

yyaxis right
sc62Y=nanmean(vmat62,1);
sc62Y=sc62Y./sc62Y(1); % GB added 8/30/21 to check normalization
drawShadedErrorRegion(binsInuM62,sc62Y,errvals62,colors{1}(2,:));
plot(binsInuM62,sc62Y,'-','Color',colors{1}(2,:), 'LineWidth',2);
ylabel('CFP/YFP Sensor Mean FRET Ratio');
ylim([.6 1]);

legend({'TomKat Sensor','CFP/YFP sensor'})
xlim([0 18])
title('Mean Cdc42 Activity v Distance from Protrusion', 'FontSize', 16)
xlabel('Distance from Protrusion (\muM)')
%%
fileNameMat=sprintf('Combo FrontBack Plot-20191121.fig');
fileNamePDF=sprintf('Comb FrontBack Plot-20191121.PDF');
savefig(g,[saveRoot filesep fileNameMat]);
printFigurePDF([saveRoot filesep fileNamePDF],'-bestfit');

%% export data in excel sheet
% loaded: SC147 FrontBack graph data.mat and SC62 FrontBack graph data.mat saved 4/4/20
%note binsInuM62=binsInuM147 will only report one.


%truncate the data at 18um.
cutInd=(binsInuM147<=19);
binsInuM2Save=binsInuM147(cutInd);
errvalsTK=errvals147(cutInd);
errvalsCY=errvals62(cutInd);
vMatTK=vmat147(:,cutInd);
vMatCY=vmat62(:,cutInd);
vMatTKmean=mean(vMatTK,'omitnan')';
vMatCYmean=mean(vMatCY,'omitnan')';


%rotate for building a table
binsInuM2Save=binsInuM2Save';
errvalsTK=errvalsTK';
errvalsCY=errvalsCY';
vMatTK=vMatTK';
vMatCY=vMatCY';

build a tk table
tTK=table(binsInuM2Save,errvalsTK,vMatTKmean,'VariableNames',{'x-Axis data (microns)', 's.e.m.','meanCellData'});
ColLabelsTK=cell(1,size(vMatTK,2));
for i=1:length(ColLabelsTK)
    ColLabelsTK{i}=sprintf('CellTK%i',i);
end
tTK1=array2table(vMatTK,'VariableNames',ColLabelsTK);
tTKFinal=[tTK tTK1];
saveLocation='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\SC147-FrontBack\Analyzed Data\2020-04-04_FB-Final-Data for both sensors';
fileNameTK=[saveLocation filesep 'frontBackDataTKSourceData.xlsx'];
writetable(tTKFinal,fileNameTK,'Sheet',1,'WriteVariableNames',true);
    
% build a cy table

tCY=table(binsInuM2Save,errvalsCY,vMatCYmean,'VariableNames',{'x-Axis data (microns)', 's.e.m.','meanCellData'});

ColLabelsCY=cell(1,size(vMatCY,2));
for i=1:length(ColLabelsCY)
    ColLabelsCY{i}=sprintf('CellCY%i',i);
end
%% plot data again to double check
figure; hold on;
yyaxis left
normTKmean=vMatTKmean./vMatTKmean(1);
plot(binsInuM2Save,normTKmean);
 drawShadedErrorRegion(binsInuM2Save,normTKmean,errvalsTK);
 ylim([.9 1]);

yyaxis right 
normCYmean=vMatCYmean./vMatCYmean(1);
plot(binsInuM2Save,normCYmean);
drawShadedErrorRegion(binsInuM2Save,normCYmean,errvalsCY);
ylim([0.6 1]);

