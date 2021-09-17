function [dat]=adjTimeStampsWithEmpiricalTimeMeasurements(dat,varargin)
% 2021-03-37 GB: The time bars for PDL plos start at t=0.3 
%for data that was collected before the metadata while the time bar made
%with metadata time stamps start at t=2. I will keep the x axis relative
%to the time of the frame before stim. Instead I will use empirical
%measurmetns to adjust the x axis for the frap time and the frame11. Using
%the measured times: Frame10 is t=0, FrapTime is t=2, Frame11 is t=4,
%frame12 is t=5.5. Based on this the stim takes ~2 sec; the frame right
%after also takes ~2sec. Then the 1.5 frame rate resumes. To avoid
%rerunning all of the data I will adjust the pDatPDL with this stim time
%paradigm. 

%added in the two pulse exp as well

% process varargin
opt.frapdelay=2;% seconds
opt.tkdelay=2; %seconds

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% determine if the dat used metadata or not
for i=1:length(dat)
    if isfield(dat(i).time,'twoPulse')
    twoPulse=dat(i).time.twoPulse;
    else
        twoPulse=false;
    end
    if ~isfield(dat(i).time,'tkInSec') && ~isempty(dat(i).time.frapTime) %determine if exp used metaData
        if twoPulse && numel(dat(i).time.frapTime)>1
            %define frames between stims
            tkTime=dat(i).time.tkTime;
            numBFind=find(tkTime>0,1);
            indBF2ndstim=find(tkTime<dat(i).time.frapTime(2),1,'last');
            tkTime(numBFind)=opt.frapdelay+opt.tkdelay;
            dat(i).time.frapTimeAdj(1)=opt.frapdelay;
            
            %define the time between stims
            betweenStimFRange=numel(tkTime(numBFind+1:indBF2ndstim));
            tkBetweenTime=[tkTime(numBFind):1.5:tkTime(numBFind)+1.5*betweenStimFRange];
            tkTime(numBFind:indBF2ndstim)=tkBetweenTime;
            dat(i).time.frapTimeAdj(2)=tkTime(indBF2ndstim)+opt.frapdelay;
            
            %define time afterStims.
            indAfter2ndStim=indBF2ndstim+1;
            afterStimRange=numel(tkTime(indAfter2ndStim+1:end));
            tkTime(indAfter2ndStim)=tkTime(indBF2ndstim)+opt.frapdelay+opt.tkdelay;
            tkAfterTime=[tkTime(indAfter2ndStim):1.5:tkTime(indAfter2ndStim)+1.5*afterStimRange];
            tkTime(indAfter2ndStim:end)=tkAfterTime;
            dat(i).time.tkTimeAdj=tkTime;
            
        else
            tkTime=dat(i).time.tkTime;
            numBFind=find(tkTime>0,1);
            tkTime(numBFind)=opt.frapdelay+opt.tkdelay;
            frameRange=numel(tkTime(numBFind+1:end));
            tkTime(numBFind:end)=[tkTime(numBFind):1.5:tkTime(numBFind)+1.5*frameRange];
            dat(i).time.tkTimeAdj=tkTime;
            dat(i).time.frapTimeAdj=opt.frapdelay;
            
        end %if twoPulse
    else
        dat(i).time.tkTimeAdj=dat(i).time.tkTime;
        dat(i).time.frapTimeAdj=dat(i).time.frapTime;
    end %if isfield(dat(i).time,'tkInSec') && ~isempty(dat(i).time.frapTime)
end %for i
end