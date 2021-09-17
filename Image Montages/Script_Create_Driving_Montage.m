%Note this script has been updated to remove nans from the stage positions
%list that are added by the multiple stimulations


%2021-04-06 I couldn't find the script that I used to make the center
%stimulation figure from fig 5a,b. I think it was this one. Sean asked me
%to update fig 5b to be relative change in fret ratio rather than the
%difference in fret ratio. I will give it a shot here.
%% set roots
masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\SC147-FrontBack\';
%masterRoot='C:\Users\George\Documents\temp Data Folder\';
  % This should be the path to the folder containing all the data from the experiment

%% camera and donut corrections and FFC
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat','donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');
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

%% imaging notes. 
% The microscope experimental script was written where the fist numbefore
% image was centered on the middle of hte cell when the image was
% collected. That image was then fed into the stageGOTOCelledgepoint
% function. Based on this, I think the first image should be disregarded
% and the pos1 values moved up. Perhaps the best way is to just pad the
% top of pos1 with zeros. In all cases the stageGOTOcellEdge is happening
% after imaging so the last pos1 value is also not useful.

%% Make image montage with FRET and fixed reference frame
rawDataDir=getDirectories([masterRoot 'Raw Data'],'2017-11-16-PLB-PP1TK_FB-older data' );
%rawDataDir=getDirectories([masterRoot 'Raw Data'],'2020-03-05_PLB_C10_driving' )
countSubdir=0; % counter for total number of subdirectories read, allows me to build a big struct for 3 days of exp
binSize=1;
imgName1='TIRF-TomKat-1';
imgName2='TIRF-TomKat-2';
numBefore=6;
for numRawDir=1:length(rawDataDir)

    clear dataSubdir; clear keepList1;
    clear p; clear frapMask;  clear goodCells;
    parentFolder=rawDataDir{numRawDir};
    mkdir([masterRoot 'Scripts' filesep parentFolder]); % makes a folder for the exp day in the scripts folder
    mkdir([masterRoot 'Analyzed Data' filesep parentFolder]);% makes a folder for the exp day in the Analyzed Data folder
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    dataSubdir=getSubdirectories(rawDataroot,'Cell31'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir);
    %load([scriptRoot filesep 'keepList1.mat'], 'keepList1');
    load([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
    load([scriptRoot filesep 'frapMask1PXrt.mat']); % now has all the normal and right shifted frapMasks in one mat file
    

    % make good cell vector
%     keepList1(5)=[];
%     goodCells=dataSubdir(keepList1);
%     
goodCells=dataSubdir;
    % iterate through the well folders
   sharpenparam=[4 25];
    mincellsize=1000;
    gRad=1;
    erodeSize=1;
    cropsize=[951 921];
    for k=1:length(goodCells)
    
        % Read stage positions
        folder=[rawDataroot filesep goodCells{k}];
        objMag=60;
        [pos1, pos1PassFlag]=createPos1(folder,objMag,rotationMatrix);
        
        % build centerMask
        % identify the cell at the center (use for frap imaging only)
            centerMask=imdilate(frapMask1PXrt, strel('disk',100));
        
        % Read and process images
        filenames1=getFilenames(folder,imgName1);
        filenames2=getFilenames(folder,imgName2);
        
        %filenames1(1)=[];
        %filenames2(1)=[];
        
        clear ratioImage; 
        %allocate ratioIm and ims cell array
        ratioImage=cell(1,length(filenames1));
        correctedRatioIm=cell(1,length(filenames1));
        grayImages=cell(1,length(filenames1));
          
        for i=1:length(filenames1)
            
            im1=double(imread([folder filesep filenames1{i}]));
            im1=(im1- medianDark1).*donutCorrection1.*halfCorrection1;
            im2=double(imread([folder filesep filenames2{i}]));
            im2=(im2-medianDark2).*donutCorrection2.*halfCorrection2;
            im = singleIms2stack(im1,im2);
            regIm=registerImagesFromQuadFit(im,p);
            cropIm2=cropImMidOut(regIm,'cropsize', cropsize);
            
            % background subtraction
            % older sean strategy
            [maskFinal, cRatioIm, imGray, acceptorIm, donorIm,minAreaThreshPass]= processFretRatioImTIRF(cropIm2,centerMask,ratioCorrectionFinal,'sharpparam',[4 50],'usemultithresh',false);

            % ratio image
            mFinal{i}=maskFinal;
            correctedRatioIm{i}=cRatioIm;
            grayImages{i}=imGray;
        end
        
        
        % determine min and max fret range for each cell
        minFret=cellfun(@(x) prctile(vect(x),1),correctedRatioIm);
        maxFret=cellfun(@(x) prctile(vect(x),99),correctedRatioIm);
        fretRange = round([min(minFret)  max(maxFret)],2);
        
        % Determine gray scale range
        
        maxGray=cellfun(@(x) prctile(vect(x),98),grayImages);
        minGray=cellfun(@(x) prctile(vect(x),5),grayImages);
        grayRange=round([mean(minGray) mean(maxGray)]);
        % Make movie frames
        
        foldNameFret=sprintf('Fret-%s',goodCells{k});
        montagePathFret=[saveRoot filesep foldNameFret];
        foldNameGray=sprintf('Gray-%s',goodCells{k});
        montagePathGray=[saveRoot filesep foldNameGray];
        mkdir(montagePathFret)
        mkdir(montagePathGray)
        
        %need to go to bed, crop the image by usign the bounding box for
        %the first and last images in the set. 
        
        % make movie for cell31
        imSelector=[1:35]; % for movie
        [xRange, yRange, ~]=getXYRangeForImMontageCrop(correctedRatioIm,imSelector,...
            'padTopLeft',40,'padBotRight',40,'pos1',pos1);
     
        
        moviePath=[saveRoot filesep goodCells{k} '.avi'];
         fretRange =[.8 1.03];
        
        writeFRETMoviesTrackedStimFromImageSequences(correctedRatioIm(1:35),...
            grayImages(1:35), moviePath,pos1,'frapspot', frapInd,'imframerate',5,...
            'numbefore',numBefore,'objective',60,'boundsfret',fretRange,...
            'xrange',xRange,'yrange',yRange,'fcounter',1);
        
     
        %make montage for cell31
        %first definve folder to save to.
        foldNameFretComposite=sprintf('Fret-%s',goodCells{k});
        montagePathFretComposite=[saveRoot filesep foldNameFretComposite];
        mkdir(montagePathFretComposite)
        % get the crop range
         imSelector=[6 8 12, 22]; % for fig 1e frames.
        [xRange, yRange, ~]=getXYRangeForImMontageCrop(correctedRatioIm,imSelector,...
            'padTopLeft',15,'padBotRight',15,'pos1',pos1);
        
        frapSelect=[6,8,12,22];
        buildCompositeFRETImgFigure(correctedRatioIm,...
            montagePathFretComposite,frapSelect,'pos1',pos1,'frapspot',frapInd,'flip',false,'numbefore',numBefore,'objective',60,...
            'xrange',xRange,'yrange',yRange,'boundsfret',fretRange,'savefiletype','pdf','imsgray',grayImages,'imframerate',5);
        
        
        
        
        
        
        
        fprintf('\nDone with set #%i.\n',k);
    end % k goodcells
end %numRawDir

%% check if frap spot is being scattered in the correct spot.
%frapcoors=frapInd+pos1 + [1-xRange(1) 1-yRange(1)];
frapcoors=frapInd+pos1;
sampleIm=correctedRatioIm{12};
frameNum=12;

Im=shiftMatrix(sampleIm,pos1(frameNum,1),pos1(frameNum,2));
imagesc(Im); hold on;
ind=[0 1 2];
colors={'r','y','p'};
for i=1:3
    
    scatter(frapcoors(frameNum+ind(i),1),frapcoors(frameNum+ind(i),2),colors{i},'LineWidth',2);
    
end
%% write a avi movie of cell driving.
bounds=[360 640 300 750] - [199 199 299 299];
%cRange=[1.8 2.05]; Cdc42
cRange=[1.9 2.4]; %old range was [1.6 1.95]

for k=1:length(goodCells)
    % * Experimental cells
    figure; axis image; set(gca,'XTick',[],'YTick',[]);
    vw=VideoWriter([root filesep goodCells{k} ' fixed frame movie_time.avi']);
    vw.FrameRate=8;
    vw.set('Quality',100);
    open(vw);
    % Read stage positions
    folder=[rawDataroot filesep goodCells{k}];
    fid=fopen([folder filesep 'stage_positions.txt'],'r');
    lines=readAllLines(fid);
    fclose(fid);
    pos1=nan(length(lines),2);
    pos1(:,1)=str2double(getTokens(lines,'x\=([\-0-9\.]+)\,'));
    pos1(:,1)=pos1(:,1)-pos1(2,1);
    pos1(:,2)=-1*str2double(getTokens(lines,'y\=([\-0-9\.]+)'));
    pos1(:,2)=pos1(:,2)-pos1(2,2);
    pos0=pos1;
    % Convert stage shifts in microns to image shifts in pixels
    pos1=pos0/distanceScale;
    pos1=(rotationMatrix*pos1')';
    pos1=round(pos1);

    % Read and process images
    filenames1=getFilenames(folder,imgName1);
    filenames2=getFilenames(folder,imgName2);
    clear images; clear ratioImage; clear imC1; clear imR1; clear im;
   
    for i=1:length(filenames1)
        
        % load images and run appropriate camera corrections
        clear im;
        im(:,:,1)=double(imread([folder filesep filenames1{i}]));
        im(:,:,1)=(im(:,:,1)-medianDark1).*donutCorrection1.*halfCorrection1;
        im(:,:,2)=double(imread([folder filesep filenames2{i}]));
        im(:,:,2)=(im(:,:,2)-medianDark2).*donutCorrection2.*halfCorrection2;
        im=registerImagesFromQuadFit(im,p);
        im=im(200:800,300:900,:);
        images{i}=im;
        
        % background subtraction
        imSum=sum(im,3);
        [maskBGinitial,objectRectangles]=fretGetInitialObjectsAndBGmask(imSum);
        imForSeg0=imageSubtractBackgroundWithObjectMaskingAndPadding(imSum,100,50,maskBGinitial);
        [maskFinal,coors]=fretGetCellMasks(imForSeg0,objectRectangles,100);
        CFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,2),100,50,maskBGinitial);
        YFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,1),100,50,maskBGinitial);
        
        % background subtraction
                    % older sean strategy
                    imSum=sum(im,3);
                    [maskBGinitial,objectRectangles]=fretGetInitialObjectsAndBGmask(imSum); % To get accurate boxes, you need to use the updated function
                    imForSeg0=imageSubtractBackgroundWithObjectMaskingAndPadding(imSum,100,50,maskBGinitial);
                    [maskFinal,coors]=fretGetCellMasks_63x(imForSeg0,objectRectangles,100);
                CFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,2),100,50,maskBGinitial);
                YFPFinal=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,1),100,50,maskBGinitial);
        
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
        ratioImage{i}=YFPsmooth./CFPsmooth;
        ratioImage{i}=ratioImage{i}./ratioCorrection(200:800,300:900);
    end
    % Register the images
    clear imC1;
    for i=1:length(filenames1)
        imC1{i}=shiftMatrix(sum(images{i},3),pos1(i,1),pos1(i,2));
        imR1{i}=shiftMatrix(ratioImage{i},pos1(i,1),pos1(i,2));
        imR1{i}(imR1{i}==0)=nan;
    end
    % Make movie frames
    for i=1:length(imC1)
        imGray=imC1{i}(bounds(1):bounds(2),bounds(3):bounds(4));
        imGray=(imGray-100)/(prctile(imGray(:),99.9)-100);
        imGray=max(0,min(1,imGray));
        imGray=repmat(imGray,[1 1 3]);
        
        imRatio=imR1{i}(bounds(1):bounds(2),bounds(3):bounds(4));
        imRatio(imRatio<cRange(1)+.01)=cRange(1)+.01;
        imRatio(isnan(imRatio))=cRange(1);
        imRatio1=grs2rgb(mat2gray(imRatio,cRange),parula_black(256));
        imDisp=  [imGray ones(size(imGray,1),2,3) imRatio1]; %add oommented region to make grayscale and fret movies
        imshow(imDisp);
        hold on
         rectangle('Position',[750 225 numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
           'LineWidth',3);
        text(780,255,sprintf('25 \x3bcm'),'Color',[1 1 1],'FontSize',16);
        hold off
        if i>6 & i<length(imC1)
            hold on;
%            scatter([frapCenter(1) + 1 -bounds(3) + pos1(i+1,1),frapCenter(2)-bounds(1)+1+pos1(i+1,2),75,[1 0.5 0],'LineWidth',2);
            scatter(frapCenter(1) + 1 -bounds(3) + [pos1(i+1,1) pos1(i+1,1)+size(imRatio,2)+2],frapCenter(2)-bounds(1)+1+[pos1(i+1,2) pos1(i+1,2)],75,[1 0.5 0],'LineWidth',2);
            hold on;
            scatter(frapCenterUsed(1) + 1 -bounds(3) + [pos1(i+1,1) pos1(i+1,1)+size(imRatio,2)+2],frapCenterUsed(2)-bounds(1)+1+[pos1(i+1,2) pos1(i+1,2)],75,[1 0 1],'LineWidth',2);
        end
        text(10,18,sprintf('t=%i sec',2*(i-6)),'Color',[1 1 0],'FontSize',14);
               
            
       
        hold off
        writeVideo(vw,getframe(gca));
        pause(0.1);
    end
    close(vw);
    close all;
    fprintf('\nDone with set #%i.\n',k);
end

%% Export colorbar

colormap('parula_black');
cbar_handle = colorbar;
printFigureEPS([root filesep 'colorbar.eps']);
imshow(imDisp) 