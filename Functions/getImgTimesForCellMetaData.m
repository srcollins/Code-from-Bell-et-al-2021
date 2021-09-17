function ts = getImgTimesForCellMetaData(metaFrame, filenames1,numBefore,varargin)
%2020-07-06 GB: 
%this function uses the compiled metadata structure (metaFrame) and pulls
%out the timeAfter data for each frame for a single cell. This data is
%stored in the structure ts. ts contains the time sequence for frames in
%the original nowInSec format (ts.tkInSec) and the elapsed time format
%(ts.tkTime) where the time from the numbefore frame was subtracted from
%all of the nowInSec time vals in the data set. Additionally, the function
%identifies the location of the frap image using the file name range
%specified by the user. The time for the frap image is also stored in TS in
%both formats.


%2020-07-30 updates: I ran into an issue where the metaData was incomplete
%for one of my experiments. I decided to add a system that checks if the
%correct frame ID was detected. If not, then the function will trigger a
%fail flag and return. I can then use the failFlag in the processing script
%to determine if the cell should be skipped.

%2021-02-20 updates: The metaData saving function on the microscope is
%currently naming images that were captured with the sequential pip3 cube
%with the same frame number, but creating two events in the meta data log.
%This is not the case for fret images that were captured simultaneously on
%two different cameras. I will add a user option for these sequential
%rfp/irfp images. Additionally, the script will check to make sure that the
%two frames are sequential in the metaFrame.

%2021-02-24 Updates: I want to use this function for experiments where
%there was no opsin stimulation. Here I will use the convention that
%numbefore=0 if there was no stim. The updated function will check to make
%sure that frap stim is >2.

%2021-03-20 Update required: Currently the function subtracts the numbefore
%from the fret image time from all fret im times. This sets the stimulation
%time to a positve number. Instead, the stim time should be subracted from
%the image times so that the stim time is 0. Im fixing in the code for now.
%% process varargin
opt.sequentialimgtf=false;

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% meat and potatoes

% pull out frame identifier from file name
fId=regexp(filenames1,'(?<=img_)[\w]{6}','match');
fId=cellfun(@(x) str2double(x),fId);

fInds=zeros(length(fId),1);
for f=1:length(fId)
    tempInd=find(metaFrame.frameID==fId(f));
    if isempty(tempInd)
        ts.failFlag(f)=true;
        ts.tkTime=inf(1,length(fInds));
        return
    else
        
        if opt.sequentialimgtf && numel(tempInd)==2
            if diff(tempInd)==1
                fInds(f)=tempInd(1);
                ts.failFlag(f)=false;
            else
                ts.failFlag(f)=true;
                ts.tkTime=inf(1,length(fInds));
                return
            end
        else
            fInds(f)=tempInd;
            ts.failFlag(f)=false;
        end
    end
end

% find frap frame to extract frap time
frapInd=find(boolRegExp(metaFrame.channel(fInds(1):fInds(end)),'FRAP'));
frapInd=frapInd+fInds(1)-1;
ts.frapInSec=metaFrame.timeAfter(frapInd);

% get tk image series timeAfter values and Id the time value for the
% numBefore frame
tAfter=metaFrame.timeAfter(fInds);
ts.tkInSec=tAfter;
if numBefore>2
    tNumBefore=tAfter(numBefore);
    tAfterElapsed=tAfter-tNumBefore;
    ts.tkTime=tAfterElapsed;
    ts.frapTime=metaFrame.timeAfter(frapInd)-tNumBefore;
    ts.timesInfrapRef=tAfter-ts.frapInSec;
else
    tNumBefore=tAfter(1);
    tAfterElapsed=tAfter-tNumBefore;
    ts.tkTime=tAfterElapsed;
    ts.frapTime=metaFrame.timeAfter(frapInd)-tNumBefore;
    ts.timesInfrapRef=[];
end
%subtract the numBefore frame from all frames to get the elapsed time
%stamps
% tAfterElapsed=tAfter-tNumBefore;
% ts.tkTime=tAfterElapsed;
% ts.frapTime=metaFrame.timeAfter(frapInd)-tNumBefore;
% ts.timesInfrapRef=tAfter-ts.frapInSec;



end
