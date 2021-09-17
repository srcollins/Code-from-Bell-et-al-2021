%% TomKat Centroid Data Graphing Script
% This script is modified from the Find the Dip script.
%2020-08-04 GB: 
%1)update res to reflect the new format.
% 2)update truncation metrics based on cell centroid movement and andmask area
% check responder v non-responder filter
% plot chuncky bins

% Used 2021-01-20 to load res20201008v3. 
% 2021-02-19 Notes: This script is still quite messy. I can be broken into
% 3 main sections. 1) Combining and cleaning up res from multiple exp. Much
% of the challenge in this section is quickly streamling the data in the
% ratioDiff sections. Currently we use the chunkybin analysis making the
% ratioDiff portions extraneous and time consuming. 2) the section uses res
% to load each imDat to apply the bleaching correction. These files are
% then used for computing the chunkybin data and the plot Data stored in the variable 
%pDat. This section of the script is the most useful/relevant. 3) plotting
%the data stored in pDat. THis section is redundant to 2021-01-23 Script to
%make centerstim plots from res and pdat. I will slim this down to make it
%more user friendly.

%2021-03-26 Modified to make photobleach supplemental figure for 15mw stim
%power. The photobleach block uses time stamps from the metadata.
%Unfortunately, the 15mw data is older and does not have metadata. 

%Note see 2021-09-15_Script to make Centerstim Plots from Res and pDat_v1.m
%to make plots from this data.

%This script will calculate a data structure called res that I use to keep
%track of each indiviual cell for all experiments. Additionally, the script
%calculates a data structure called pDat that compiles the data from cells of
%like groups for plotting. 
%% set root
%masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim\';
masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-CenterStim-Fast\';

%masterRoot = 'C:\Users\George\Documents\temp Data Folder\';  % This should be the path to the folder containing all the data from the experiment
resRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim\';


%% build Res from saved imDats
res=createResFromImgDat(masterRoot,'dirkeyword','Fast-M');

%% save Res
tempSaveRoot='C:\Users\George\Documents\temp Data Folder\Analyzed Data\Fast res';
FastRes5Pulse20201028_v1=res;
save([tempSaveRoot filesep 'FastRes5Pulse20201028_v1.mat'],'FastRes5Pulse20201028_v1', '-v7.3');

%% load res data structure
%load([root filesep 'Compiled_Cdc42TK_Centroid_Data1PX-LatA-12pX.mat'],'res1PXrtLatA');
% load([saveRoot filesep 'centerStimRes20200625.mat'],'centerStimRes20200625');
% load([saveRoot filesep 'centerStimRes20200708.mat'],'centerStimRes20200708');
% load([saveRoot filesep 'centerStimRes20200710and11.mat'],'centerStimRes20200710and11');
% load([saveRoot filesep 'centerStimRes20200715.mat'],'centerStimRes20200715');
% load([saveRoot filesep 'centerStimRes20200716To22.mat'],'centerStimRes20200716To22');
% bigRes=[centerStimRes20200625 centerStimRes20200708 centerStimRes20200710and11 centerStimRes20200715,centerStimRes20200716To22];

%load([saveRoot filesep 'res20200902v3.mat'])
%% clean up bigRes

bigRes=cleanUpBigRes(bigRes,masterRoot);
%% Add the time data saved in imDat.ts to bigRes

for i=1:length(bigRes)
    fileNamePath=assembleFilePathFromCentStimStuct(masterRoot,bigRes,i);
    
    for j=1:length(fileNamePath)
          load(fileNamePath{j});
          frapTimeAtZero=imDat.times.tkInSec-imDat.times.frapInSec(1);
    bigRes(i).times{j}=imDat.times.tkTime';
    bigRes(i).frapTimeAtZero{j}=frapTimeAtZero;
    end
end



%% Distance Scale
umPerPix = distanceScale(60);  %=0.21 um/pixel for 60x


%% save bigRes W 12 pix movement thresh. Pre non-responder removal
ResSaveRoot=[resRoot 'Analyzed Data' filesep '2020-07-30_PooledDatafromJun2020-july2020' ];
bigRes20200825W12pixThreshAllcells=bigRes;
save([ResSaveRoot filesep 'bigRes20200825W12pixThreshAllcells.mat'],'bigRes20200825W12pixThreshAllcells','-v7.3');

%% define res
%res=res20200121Cdc42KO_PakKOAllredo;
res=bigRes;
%res=rescdc42ko2020CenterStim;
%% implement the frap stim photobleaching correction
% This section loads the bleacing correction parameters and builds a
% structure array that contains the nominal and actual frap measurements.
% The multiPlier value was empirically derived.

pBleachCorrRoot='D:\Users\collinsLab\Documents\GB_Data\Image Corrections\2020-08-28-15mw-FRAP-Photobleaching Corrections\';
load([pBleachCorrRoot 'Local Bleaching Correction Parameters 2020-08-28.mat']);
multiPlier=2;
nominalFrap={'2mw','3mw','5mw','10-5mw','15mW'}; %nominal laser poer 
actualFrap=[1.78, 4.35, 9.50, 0.8,37.58]*multiPlier;
bleacPwrNom={'15'};
bleachPwrAct=37.58;
for i=1:length(nominalFrap)
    frapPwr(i).bleacPwrNom=bleacPwrNom;
    frapPwr(i).bleachPwrAct=bleachPwrAct;
    frapPwr(i).nominal=nominalFrap{i};
    frapPwr(i).actual=actualFrap(i);
end



for r=1:length(res)
    fileNamePath=assembleFilePathFromCentStimStuct(masterRoot,res,r);
    for f=1:length(fileNamePath)
        clear imDat distMask bleachCorrectionImage;
        load(fileNamePath{f});
        %calc time stamps and find time inds that are post stim
        timeStamps=imDat.times.tkInSec-imDat.times.frapInSec(1);
        postStimInd=timeStamps>0;
        postStimInd=postStimInd'; 
        % determine nominal laserpower from file name
        nomLPFromFileName=regexp(imDat.fileName,'[0-9]+(?=mW)','match'); % may need to adjust for filenames. If all exp are 2mw just set equal to '2mW';
        nomLPFromFileName=strcat('^',nomLPFromFileName);
        nomLPind=find(arrayfun(@(x) boolRegExp(x.nominal, nomLPFromFileName{1}),frapPwr));
        % build relativeLaserPower vect
        relativeLaserPower=repmat(frapPwr(nomLPind).actual/frapPwr(i).bleachPwrAct,1,length(postStimInd));
        relativeLaserPower=relativeLaserPower.*postStimInd;
        
        if length(imDat.times.frapTime)<2
        % build distMask image
        for img=1:length(imDat.cropFret)
        %distMask{img}=bwdistgeodesic(logical(imDat.cropMask{img}),logical(imDat.cropFrapIm{img}),'quasi-euclidean')*distanceScale(60);
        distMask{img}=bwdist(logical(imDat.cropFrapIm{img}))*distanceScale(60);
        bleachCorrectionImage{img}=makeLocalBleachCorrectionImage(pBleachCorr,double(distMask{img}),relativeLaserPower(img),timeStamps(img));
%         imDat.pBleachCropFret{img}=imDat.cropFret{img}./bleachCorrectionImage{img};
        imDat.bleachCorrection{img}=bleachCorrectionImage{img};
        %imDat.postStimInd=postStimInd;
        end
        save(fileNamePath{f},'imDat','-v7.3');
        else
            for img=1:length(imDat.cropFret)
                distMask{img}=bwdist(logical(imDat.cropFrapIm{img}))*distanceScale(60);
                bleachCorrectionImage{img}=makeLocalBleachCorrectionImage(pBleachCorr,double(distMask{img}),relativeLaserPower(img),timeStamps(img));
            end
            compositeBleachCorr=buildCompositeBleachCorrectionForMultiStim(bleachCorrectionImage,imDat);
            imDat.bleachCorrection=compositeBleachCorr;
            save(fileNamePath{f},'imDat','-v7.3');
        end
    end
end



%% Run chunkyBin Analysis .
%This is the center stim analysis strategy described in Figure 7c.

fCount=0;
% fRInitial=[1:60];
% frameRange=[1:60];
clear pDat sd sem sdAll semAll yVals normY allNormY condLabel stimLabel
for f=1:length(res)
     times=[];
     timesInfrapRef=[];
    for t=1:length(res(f).times)
        times=[times;res(f).times{t}];
       timesInfrapRef=[timesInfrapRef; res(f).frapTimeAtZero{t}'];
    end
    meanTimes=round(nanmean(times),2);
    meanTimeFrapRef=round(nanmean(timesInfrapRef),2); 
    
    fileNamePath=assembleFilePathFromCentStimStuct(masterRoot,res,f);
    fCount=fCount+1;
    clear ratioDat cBinOut
    for o=1:length(fileNamePath)
        load(fileNamePath{o});
          %[cBinOut]= calcChunkyBinsV2(imDat.cropKat,imDat.cropTom,imDat.cropFrapIm,imDat.cropRatioCorrection,'binedges', [1:1:12 ]);
             [cBinOut]= calcChunkyBinsV3(imDat.cropKat,imDat.cropTom,...
                imDat.cropFrapIm,imDat.cropRatioCorrection,imDat.bleachCorrection,'binedges', [1:1:12 ]);
           
            numCurves=find(boolRegExp(fieldnames(cBinOut(1)),'ratioBin'));
            areaCheckInd=find(boolRegExp(fieldnames(cBinOut(1)),'binArea'));
            fieldNAll=fieldnames(cBinOut);
            
            
            for nc=1:length(numCurves)
                % get ratio for all frames in each bin
                fieldN=fieldNAll{numCurves(nc)};
                tempVect=arrayfun(@(x) x.(fieldN),cBinOut);
                %check area
                fNameA=fieldNAll{areaCheckInd(nc)};
                tempAreaV=arrayfun(@(x) x.(fNameA),cBinOut);
                tempAreaV(isnan(tempAreaV))=0;
%                 ratioDat.(fNameA)(o,:)=tempAreaV(fRInitial);
%                 ratioDat.(fieldN)(o,:)=tempVect(fRInitial);
                ratioDat.(fNameA)(o,:)=tempAreaV;
                ratioDat.(fieldN)(o,:)=tempVect;
                if any(isnan(ratioDat.(fieldN)(o,2:end))) || any(ratioDat.(fNameA)(o,2:end)<cBinOut(1).areaThresh)
                    ratioDat.(fieldN)(o,:)=nan;
                end
            end %for nc
        
    end  %for o
    
    [ratioDat,numRm,numCellTot]=removeCellsChunkyBinData(ratioDat,'notbin0',false);
    ratioDatFName=fieldnames(ratioDat);
    rDatInd=boolRegExp(ratioDatFName,'ratioBin');
    ratioDatFName=ratioDatFName(rDatInd);
    
    colors=jet(length(ratioDatFName));
    for j=1:length(ratioDatFName)
        
        %sd=nanstd(ratioDat.(ratioDatFName{j})(:,frameRange));
        sd=nanstd(ratioDat.(ratioDatFName{j})(:,:));
        %sem=nanstd(ratioDat.(ratioDatFName{j})(:,frameRange))./sqrt(size(ratioDat.(ratioDatFName{j})(:,frameRange),1));
        sem=nanstd(ratioDat.(ratioDatFName{j})(:,:))./sqrt(size(ratioDat.(ratioDatFName{j})(:,:),1));

        sdAll(j,:)=sd;
        semAll(j,:)=sem;
        %yVals=nanmean(ratioDat.(ratioDatFName{j})(:,frameRange));
        yVals=nanmean(ratioDat.(ratioDatFName{j})(:,:));
        normY=yVals./yVals(2);
        allNormY(j,:)=normY;
      
        
        
        condLabel=regexp(res(f).name{1},'.*(?=-SC344)','match');
        stimLabel=regexp(res(f).name{1},'(?<=centroid_)[a-zA-Z0-9]+','match');
        
        pDat(fCount).stimLabel=stimLabel; 
        pDat(fCount).condLabel=condLabel;
        pDat(fCount).label=res(f).name{1};
        pDat(fCount).sd=sdAll;
        pDat(fCount).sem=semAll;
        pDat(fCount).rawY=ratioDat.(ratioDatFName{j});
        pDat(fCount).normY=allNormY;
        %pDat(fCount).xval=meanTimes(frameRange);%[1:length(pDat(fCount).normY)];
        pDat(fCount).xval=meanTimes;
        %pDat(fCount).xvalFrap=meanTimeFrapRef(frameRange);
        pDat(fCount).xvalFrap=meanTimeFrapRef;
        pDat(fCount).numRM=numRm;
        pDat(fCount).numCellTot=numCellTot;
        pDat(fCount).removeInd=find(isnan(ratioDat.(ratioDatFName{2})(:,2))); %indicies for cells that were set to nan due to low pixel counts
        pDat(fCount).legends=ratioDatFName;
        pDat(fCount).colors=colors;
        pDat(fCount).(ratioDatFName{j})=ratioDat.(ratioDatFName{j});
        
    end
end


    
%% save res and pDat for future use
resRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-CenterStim-Fast\Analyzed Data\Fast res';
pDat20210219pBleach=pDat;
save([resRoot filesep 'pDat20210219pBleach.mat'],'pDat20210219pBleach','-v7.3');
res20210219v6=res;
save([resRoot filesep 'res20210219v6.mat'],'res20210219v6','-v7.3');


