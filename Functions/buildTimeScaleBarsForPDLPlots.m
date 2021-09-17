function timeBars = buildTimeScaleBarsForPDLPlots(pDatPDL,inds2plot, colors, varargin)
%2020-10-21 GB: This function is designed to add stimulation durtation
%scalebars to PDL assay plots. These bars should also match the color of
%the corresponding curve.

%20201217 GB: updated this funciton to produce a time bar that is overlayed
%to take up less space. Allows user to specify overlay feature in varargin.

%2021-03-20 GB update: The timebars and plots should use the first stimulus
%time as the reference time, not the numbefore frame as it exists
%currently. I added a varargin input that will allow the user to select the
%pDatPDL.t0isStimTime option that uses the time stamps where the first stim
%is the t=0 ref. The input ('imfieldname') will be a string that contains the fieldname to
%use (eg pDatPDL.t0isStimTime). To use this alternative input the user should add a second 
% varargin input called 'stimfieldname' to also give the fieldname for the
% stim times in the stim ref (pDatPDL.stimInStimRef).

%2021-06-10 GB: For the Cdc42 paper we ended up plotting time on the x-axis
%relative to the Cdc42 FRET image immediately preceding stimulation. Thus
%the t0isStimTime should be disregarded. Instead use 'meanTKtime'. 

% Note: if you need to plot timebars for two pulse experiments use:
%buildTimeScaleBarsFor2PulsePDLPlots instead.

%% process varargin 
opt.miny=[];
opt.minheight=0.001;
opt.overlay=0;
opt.imfieldname='';
opt.stimfieldname='';
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% determine the fieldnames to use
if isempty(opt.imfieldname) && isempty(opt.stimfieldname)
    imFieldName='meanTKtime';
    stimFieldName='meanStimTime';
else
    imFieldName=opt.imfieldname;
    stimFieldName=opt.stimfieldname;
end
    


%% make time bars
% determine if any of the inds2plot are 0 stim conditions.
stimInd=[];
for i=1:length(inds2plot)
    stimInd(i)=any(~isnan(pDatPDL(inds2plot(i)).(stimFieldName)));
    %stimInd(i)=any(~isnan(pDatPDL(inds2plot(i)).meanStimTime));
end

if sum(stimInd)==0
    fprintf('Non Stimulated');
    timeBars(1).x=[];
    timeBars(1).y=[];
    timeBars(1).w=[];
    timeBars(1).h=[];
    timeBars(1).colors=[];
    timeBars(1).numStim=[];
    return
end

%preallocate timebars
for j=sum(stimInd):-1:1
   timeBars(j).x=[];
   timeBars(j).y=[];
   timeBars(j).w=[];
   timeBars(j).h=[];
   timeBars(j).colors=[];
   timeBars(j).numStim=[];
end

%% determine placement
% first find min yVal
if isempty(opt.miny)
minYdata=round(min(arrayfun(@(x) min(x.meanNormRatio,[],'all'), pDatPDL)),3);
maxY=minYdata-.005;
else 
    maxY=opt.miny;
end

minH=opt.minheight;

% update if loop here 20201217
if opt.overlay
     minY=maxY;
    yRange=repmat(minY,1,sum(stimInd));
else
    minY=maxY+((minH*2)*sum(stimInd));
    yRange=linspace(maxY, minY,sum(stimInd));
    yRange=fliplr(yRange);
end

inds2PFinal=inds2plot(logical(stimInd));

%apply stimInds to colors
stimIndColors=logical(stimInd');
%stimIndColors=logical(repmat(stimIndColors,1,3));
colors(~stimIndColors,:)=[];

% fill in timeBars
count=0;
for t=1:length(inds2PFinal)
    %get numStim for ranking later
    numStim=regexp(pDatPDL(inds2PFinal(t)).numStim,'[\d]+','match');
    numStim1=str2num(char(numStim{1}));
    timeBars(t).numStim=numStim1;
    nextTKimgInd=[];
    nextTKimg=[];
    timeBars(t).x=pDatPDL(inds2PFinal(t)).(stimFieldName)(1);
    timeBars(t).y=yRange(t);
    if sum(~isnan(pDatPDL(inds2PFinal(t)).(stimFieldName)))==1
        
        nextTKimgInd=find(pDatPDL(inds2PFinal(t)).(imFieldName) > timeBars(t).x);
        nextTKimg=pDatPDL(inds2PFinal(t)).(imFieldName)(nextTKimgInd(1));
        timeBars(t).w=nextTKimg-timeBars(t).x;
    else
        endStim=pDatPDL(inds2PFinal(t)).(stimFieldName);
        endStim(isnan(endStim))=[];
        nextTKimgInd=find(pDatPDL(inds2PFinal(t)).(imFieldName) >endStim(end));
        nextTKimg=pDatPDL(inds2PFinal(t)).(imFieldName)(nextTKimgInd(1));
        timeBars(t).w=nextTKimg-timeBars(t).x;
    end
    timeBars(t).h=minH;
    timeBars(t).colors=colors(t,:);
end

end