%% Script for computing Alignment parameters, FrapMasks, and KeepLists
%Notes: This script was used for all under-agarose preparation experiments.
%The first sub-section defines the alignment parameters saved as a data
%structure called 'p'. 

% The second sub-section defines a logical mask image of the frap spot. This mask
% is a single true pixel for the center of the frapspot, which is defined as the 
% brightest pixel in the collected frap images.

%Finally, the third subsection allows the user to visualize every
%identified cell and select single cells for later processing. The output
%is a saved variable called keepList1.

%% set roots
masterRoot='C:\Users\George\Documents\temp Data Folder\';
%masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-HighRes-CenterStim\';
%% Is this TomKat Img Analysis?
TFtomKat=true; 

%% camera and donut corrections and FFC
%load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat','donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');
load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\ratioCorrectionFinalsTK.mat','ratioCorrectionFinal24');
ratioCorrectionFinal=ratioCorrectionFinal24;

%% compute alignment and ratio correction data for each experimental day
rawDataDir=getDirectories([masterRoot 'Raw Data'],'MultiP' )


countSubdir=0; % counter for total number of subdirectories read, allows me to build a big struct for 3 days of exp
binSize=1;
% imgName1='TIRF-TomKat-1';
% imgName2='TIRF-TomKat-2';

imgName1='TomKat-1';
imgName2='TomKat-2';

for numRawDir=1:length(rawDataDir)
    clear dataSubdir;
    parentFolder=rawDataDir{numRawDir};
   mkdir([masterRoot 'Scripts' filesep parentFolder]); % makes a folder for the exp day in the scripts folder
   mkdir([masterRoot 'Analyzed Data' filesep parentFolder]);% makes a folder for the exp day in the Analyzed Data folder
    rawDataroot = [masterRoot 'Raw Data' filesep parentFolder];  % This should be the path to the folder containing all the data from the experiment
    scriptRoot=[masterRoot 'Scripts' filesep parentFolder];
    saveRoot=[masterRoot 'Analyzed Data' filesep parentFolder];
    dataSubdir=getSubdirectories(rawDataroot,'centroid'); % get the subdirectories for each folder in rawData
    dataSubdir=natsortfiles(dataSubdir); % sort the data
    
  foldDate=str2double(regexp(parentFolder, '^[0-9]+','match'));
    if foldDate < 2020
        load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera corrections 2018-08-09.mat','donutCorrection1','donutCorrection2','halfCorrection1','halfCorrection2','medianDark1', 'medianDark2');
    else
        load('D:\Users\collinsLab\Documents\GB_Data\Image Corrections\Camera correctionsTK 2020-03-23.mat');
    end
   
    
    %-- computes alignment if not completed already
    if exist([scriptRoot filesep 'alignment parameters pX pY.mat'])~=2;
        if TFtomKat
            p0=buildP0TomKat;
        else
            p0=buildP0TomKat('x1',0,'y1',0);
        end
        [p, imgF] = collectImgAndAlignWithPQuadFit(rawDataroot, dataSubdir(1:20), {imgName1,imgName2}, 'p0', p0, 'binsize',1);
        showImagesMergeChannels(imgF(:,:,1),imgF(:,:,2));
        save([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
    else
        % if P is defined, load p parameters.
        load([scriptRoot filesep 'alignment parameters pX pY.mat'],'p');
    end % if exist([scriptRoot filesep 'alignment parameters pX pY.mat'])==0;
    
    
   
    % make a frap mask. There are a few options depending on the
    % experiment. If a bright frap image was collected, use
    % make1PXFrapMask2CameraTemp to make the frap mask. Alternatively, if
    % no bright image was collected, use makeCompositeFrapIm to sum all
    % regular frap images to generate a single image that highlights the
    % frap spot position. This can be tricky because the cameras are not
    % synchronized with the frap stimulation. Finally if you know the frap
    % coordiantes, you can make a mask using
    % makeRegisteredCroppedFrapImFromCoors. 
    
    
    if exist([scriptRoot filesep 'frapMask1PXrt.mat'])~=2;
        cropsize=[951 921];  %cropRange=[200,800,300,900];
        deltaSpot = defineDeltaSpotFromStagePos(rawDataroot,dataSubdir,{imgName1,imgName2},'numbefore',10);
        
        prompt= 'Do you have a bright frap image? use true/false: ';
        resp1= input(prompt);
        if resp1
            keyWord='FRAP'; imId='FRAP'; 
             
        [frapMask1PX,frapInd,frapMask1PXrt frapIndrt]=...
            make1PXFrapMask2CameraTemp(rawDataroot,p,'keyword',keyWord,'imid',imId,'cropsize',cropsize,'deltaspot',deltaSpot);
        save([scriptRoot filesep 'frapMask1PXrt.mat'],'frapMask1PXrt','frapIndrt','frapMask1PX','frapInd');
        else % if resp1
        prompt = 'enter the frap coordinates in [x y] vect. If no coors enter: n';
        resp2=input(prompt);
          if isvector(resp2)
              [frapMask1PX,frapInd,frapMask1PXrt, frapIndrt]...
                  = makeRegisteredCroppedFrapImFromCoors(resp2, p,'deltaspot',deltaSpot,'cropsize',cropsize);
              save([scriptRoot filesep 'frapMask1PXrt.mat'],'frapMask1PXrt','frapIndrt','frapMask1PX','frapInd');
          else
              keyWord='centroid'; imId='FRAP'; cropsize=[951 921]; %cropRange=[200,800,300,900];
              sumFrapIm=makeCompositeFrapIm(rawDataroot,p,'keyword',keyWord,'imid',imId,'cropsize',cropsize); %THIS OPTION IS THE WORST ONE!!
              [frapMask1PX, frapInd] = make1PXFrapMaskWithImgInMemory(sumFrapIm,'cropsize',cropsize);
          end %isvector(resp2)
        end % if resp1
         
       
    else
        load([scriptRoot filesep 'frapMask1PXrt.mat'],'frapMask1PXrt','frapIndrt','frapMask1PX','frapInd');
    end % if exist([scriptRoot filesep 'frapMask1PXrt.mat'])~=2;
    
    
     % generate keepList1 or load a keepList if one exists
    if exist([scriptRoot filesep 'keepList1.mat'])~=2;
        cropsize=[951 921];
        keyword='centroid';
        %manual keep list
        [keepList1,keepListFoldNames]=generateKeepList(rawDataroot,dataSubdir,{imgName1,imgName2},p,keyword, 'frapcoors',frapIndrt,'cropsize',cropsize);
        save([scriptRoot filesep 'keepList1.mat'],'keepList1','keepListFoldNames','oStats');
    else
        load([scriptRoot filesep 'keepList1.mat'],'keepList1');
    end
    

    
end% numRawDir   
 
    
 %% play a quick movie to check a keepList
 {imgName1,imgName2}
 goodCells=dataSubdir(keepList1);
 for siteCheck=1:length(goodCells)
    
    clear im
    folder=[rawDataroot filesep goodCells{siteCheck}];
    numBefore = findNumImgBeforeFrapStim(folder, {imgName1,imgName2});
        filenames1=getFilenames(folder,imgName1);
        filenames2=getFilenames(folder,imgName2);
        
       
            figure(1);
            for frameNum=1:length(filenames1)
            im(:,:,1)=double(imread([folder filesep filenames1{frameNum}]));
            im(:,:,2)=double(imread([folder filesep filenames2{frameNum}]));
            im2=imageBin(registerImagesFromQuadFit(im,p),1);
            im2=cropImMidOut(im2,'cropsize', cropsize);
            showImagesMergeChannels(im2(:,:,1),im2(:,:,2));
            cellFolder=sprintf('%s',goodCells{siteCheck});
            text(20,20,cellFolder,'Color','w');
            drawnow;
            if ~isempty(frapInd) && frameNum> numBefore
                hold on; scatter(frapInd(1), frapInd(2), 'r', 'LineWidth', 3);
%                 scatter(opt.frapcoors(2), opt.frapcoors(1), 'w', 'LineWidth', 3);
            end
            pause(0.05);
            end %img play movie loop
            pause(0.1);
            
            
            close(1);
            fprintf('wells remaining:%i\n',length(goodCells)-siteCheck);
     
end %sitecheck

