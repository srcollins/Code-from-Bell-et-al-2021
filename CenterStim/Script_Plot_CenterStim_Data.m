%% 2021-01-23 Script to create centerstim plots from saved res and pDat files
%My centerstim analysis script was getting really crowded so I decided to
%make a new script that loads the saved res and pDat structures and is then
%used to generate plots for my figures.

%2021-03-25: updated to make the chuncky bin plot for SFig5b. 

%% roots
masterRoot='D:\Users\collinsLab\Documents\GB_Data\Parapinopsin Cdc42TK experiments\PP1-Centerstim-Fast\';
%% load saved res and pDat data
% drag and drop is easiest. Used the following files:
%2021-01-29 I recomputed the cell masks to be tighter using the bounding  box. This
%new data is from the Rerun of the 1 pulse fast low power data. 
%res20210129v5;
%pDat20210129pBleach;
%% rename res and pDat for ease of use
 res=res20210129v5;
 pDat=pDat20210129pBleach;

%% save a copy of this script in the scriptRoot

Version='SFig5';
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

%% Generate relative fret ratio v distance plots from chuncky bin data
%For this data set, I need to go into the ratio bins and calculate the
%relative fret (chunkyBinDataForAllFrames./chunkyBinDataForFrameNumBefore)
%then I need to calculate the sd and SEM so that I can make these plots. I
%will additional fields in pDat as part of this process so that I can have
%the data for later.

% get ratioBin field names from the pDat structure array
fieldNAll=fieldnames(pDat);
fieldInd=boolRegExp(fieldNAll,'ratioBin');
fieldNAll=fieldNAll(fieldInd);
% determine the number of frames before stim. Use pDat(1).xvalFrap and find
% # of negative vals.
numBefore=sum(pDat(1).xvalFrap<0);

% this for loop series will pull out the 
for p=1:length(pDat) % iterate through the structures in pDat array.
    relSd=nan(length(fieldNAll), size(pDat(1).ratioBin0,2));
    relSem=nan(length(fieldNAll), size(pDat(1).ratioBin0,2));
    relRatioBinFret=nan(length(fieldNAll), size(pDat(1).ratioBin0,2));
    rBinsTemp=[];
     numCells=[];
    
    for rbins=1:length(fieldNAll) % iterate through the ratioBin field names.
        rBinsTemp=pDat(p).(fieldNAll{rbins}); % pull out data for each ratioBin
        rBinsTemp=rmmissing(rBinsTemp); %remove nans
        numCells(rbins)=size(rBinsTemp,1);
        for cells=1:size(rBinsTemp,1)
            rBinsTemp(cells,:)=rBinsTemp(cells,:)./rBinsTemp(cells,numBefore); %compute fret vals relative to numbefore frame
        end
        % calculate nanmean, sd and sem for all cells in the ratioBin. It
        % is still iterating over the length of rBins so the final should
        % be a 13 x 30 matrix where the rows are the ratiobins and the
        % columns are the frames.
        relSd(rbins,:)=std(rBinsTemp);
        relSem(rbins,:)=std(rBinsTemp)./sqrt(size(rBinsTemp,1));
        relRatioBinFret(rbins,:)=nanmean(rBinsTemp);
    end
   
    pDat(p).numCellsNansRM=mode(numCells);
    pDat(p).relSd=relSd;
    pDat(p).relSem=relSem;
    pDat(p).relRatioBinFret=relRatioBinFret;
end
%% remove outer cell array for stim and cond labels if necessary
% stim and cond labels were a nested cell array for the fast 1p data. 

if iscell(pDat(1).condLabel{1})
    for p=1:length(pDat);
        pDat(p).stimLabel=[pDat(p).stimLabel{:}];
        pDat(p).condLabel=[pDat(p).condLabel{:}];
    end
end
%% plot relative fret ratio as a function of distance from stim site
% This block will plot relative fret ratios v distance using ther error bar
% scatter plot function. Its set up to plot a single pDat Condition.

%Figure 7d,g,i
%id the pDat index to plot
stimCond={'all'};
condType={'^Ctrl'};
numStim=[];
pInd=sortpDatCenterStimbyStimConds(pDat,condType,stimCond,numStim);


%determine the frame before stim
numBefore=sum(pDat(1).xvalFrap<0);

binRange=[1:10];
distXvals=[0:9]; 
% set the range for frames to be plotted
frameRange=[numBefore:2:numBefore+11];
legends={};
for j=pInd%:length(pDat)
     relFretYVals=[]; relFretSem=[];
    %need to remove bin0 as its the whole cell mean.
    relFretYVals=pDat(j).relRatioBinFret(2:end,frameRange);
    relFretYVals=relFretYVals(binRange,:);
    relFretSem=pDat(j).relSem(2:end,frameRange);
    relFretSem=relFretSem(binRange,:);
    
    colors=parula_black(length(frameRange));
    
    for i=1:length(frameRange)
        hold on;
        errorbar(distXvals,relFretYVals(:,i),relFretSem(:,i),'LineWidth',1.25,'Color',colors(i,:));
        %plot(distXvals,relFretYVals(:,i),'-o','LineWidth',1.25,'Color',colors(i,:));
        lName={};
        if i==1
            lName=sprintf('t= Before');
        else
            lName=sprintf('t= %s sec',num2str(round(pDat(j).xvalFrap(frameRange(i)),1)));
        end
        legends{i}=lName;
    end
    ylabel('Change in Signaling (Relative Fret Ratio)')
    xlabel('Dist From Target (\mum)')
    
    %ylim([0.995 max(relFretYVals,[],'all')])
    ylim([.997 1.016])
    yticks([1:0.005:1.015])
    xlim([-0.5 9.5])
    %title(sprintf('%s-%s %s Cells',pDat(j).stimLabel{1},pDat(j).condLabel{1}, num2str(pDat(j).numCellsNansRM)),'Interpreter','none')
end
legend(legends,'Location','northeast');
%% plot relative Fret v dist for specific timepoints across the three conditions.
%Figure 7e&f

stimCond={'all'};
condType={'all'};
numStim=[];
inds2plot=sortpDatCenterStimbyStimConds(pDat,condType,stimCond,numStim);

%determine the frame before stim
numBefore=sum(pDat(1).xvalFrap<0);

binRange=[1:10];
distXvals=[0:9]; 
% set the range for frames to be plotted
frameRange=[numBefore+2 numBefore+8];
colors=parula_black(length(inds2plot));

for i=1:length(frameRange)
     %subplot(1,2,i);
     figure;
    hold on; legends={};
    
    for p=1:length(inds2plot)
        % build tempmatrix to contain the relative Fret data for a specific
        % frame.
        relFretYVals=[]; relFretSem=[];
        %need to remove bin0 as its the whole cell mean.
        relFretYVals=pDat(inds2plot(p)).relRatioBinFret(2:end,frameRange(i));
        relFretYVals=relFretYVals(binRange,:);
        relFretSem=pDat(inds2plot(p)).relSem(2:end,frameRange(i));
        relFretSem=relFretSem(binRange,:);
        %plot
        errorbar(distXvals,relFretYVals,relFretSem,'Color',colors(p,:));
        
        %build legends
        lName={};
        lName=sprintf('%s',pDat(inds2plot(p)).condLabel{1});

        legends{p}=lName;
    end
    legend(legends);
   yline(1, '--', 'Color', [17 17 17]/255);
    ylabel('Change in Signaling (Relative Fret Ratio)')
    xlabel('Dist From Target (\mum)')
    
    ylim([0.997 1.022])
    yticks([1:0.01:1.02]);
    xlim([-0.5 9.5])
    title(sprintf('Spatial Response at t= %s sec',num2str(round(pDat(j).xvalFrap(frameRange(i)),1))),'Interpreter','none')
end
% legend(legends,'Location','northeast');
    



%% plot chunkybin data with fewer curves
% Creating chunky bin supplemental figures that will be for the high power
% (5mw) stim condition for Ctrl, LatA and KO conditions.

%%   subplot means for each 5mw bin for all conds(ctrl, LatA, C10)
%2021-03-25 altered to plot bin1 for all three conditions. 
%SFig 5d&e
stimCond={'10'};
condType={'^Ctrl$','^C10-Cdc','^LatA$'};
numStim=[];
inds2plot=sortpDatCenterStimbyStimConds(pDat,condType,stimCond,numStim);
pInd=inds2plot;
pRange=[5]; % sets the ratio bins to be plotted
% tempColors=pDat(1).colors;
% tempColors(9,:)=tempColors(10,:);
tempColors=parula(4);
figure; hold on;
clear legendSelect
 for i=1:3%:length(pInd)
     %sp=subplot(1,length(pInd),i);
     
     
     %build tempNormY that contains the y data
     tempNormY=[];
     tempNormY=pDat(pInd(i)).normY(2:end,:);
     %build tempSEM that contains the sem data
     tempSEM=[];
     tempSEM=pDat(pInd(i)).sem(2:end,:);
     % build temp legends to remove bin0
%      tempLegends=pDat(pInd(i)).legends(2:end);
     
     legendSelect{i}=sprintf('%s',pDat(i).condLabel{1});
     for j=1:length(pRange)
         drawShadedErrorRegion(pDat(i).xvalFrap,tempNormY(pRange(j),:),tempSEM(pRange(j),:),tempColors((i),:));
         plot(pDat(i).xvalFrap,tempNormY(pRange(j),:),'Color',tempColors((i),:),'LineWidth',1.5);
         tempLegend={};
         tempLegend=regexp(tempLegends{pRange(j)},'(?<=ratioBin).*','match');
     end
     ylabel('Cdc42 Activity (Mean Fret Ratio)','FontSize',8);
     xlabel('Time (sec)','FontSize',8);
     %title(sprintf('%s-%s %s um Bin %s Cells',pDat(i).stimLabel{1},pDat(i).condLabel{1},num2str(cBinOut(1).binWidth), num2str(pDat(i).numCellTot)),'Interpreter','none');
     title('0.8 microWatts'); 
     
     %title(sprintf('%s-%s %s Cells',pDat(pInd(i)).stimLabel{1},pDat(pInd(i)).condLabel{1}, num2str(pDat(pInd(i)).numCellTot)),'Interpreter','none');
     legend(legendSelect,'FontSize',6);
     ylim([0.99 1.035]);
     xlim([min(round(pDat(1).xvalFrap,1)), 20]);
     %xlim([min(round(pDat(1).xvalFrap,1)) max(round(pDat(1).xvalFrap,1))]);
 end
