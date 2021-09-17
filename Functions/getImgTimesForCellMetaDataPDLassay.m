function ts = getImgTimesForCellMetaDataPDLassay(metaFrame, filenames1,numBefore,varargin)
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

%2020-08-06 Updates: need to handle when numbefore is 0,1,or multiples.

%2021-03-20 Update required: Currently the function subtracts the numbefore
%from the fret image time from all fret im times. This sets the stimulation
%time to a positve number. Instead, the stim time should be subracted from
%the image times so that the stim time is 0. Im fixing in the code for now.
%% process varargin
opt.stimterm='DAPI';
opt.numbefore=10;
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
fInds(f)=tempInd;
ts.failFlag(f)=false;
    end
end

% find frap frame to extract frap time
if ~isempty(numBefore)
frapInd=find(boolRegExp(metaFrame.channel(fInds(1):fInds(end)),opt.stimterm));
frapInd=frapInd+fInds(1)-1;
ts.frapInSec=metaFrame.timeAfter(frapInd);

% get tk image series timeAfter values and Id the time value for the
% numBefore frame
tAfter=metaFrame.timeAfter(fInds);
ts.tkInSec=tAfter;
tNumBefore=tAfter(numBefore);   
ts.stimInSec=metaFrame.timeAfter(frapInd);
%subtract the numBefore frame from all frames to get the elapsed time
%stamps
tAfterElapsed=tAfter-tNumBefore;
ts.tkTime=tAfterElapsed;
ts.frapTime=metaFrame.timeAfter(frapInd)-tNumBefore;
else
    numBefore=opt.numbefore;
    frapInd=find(boolRegExp(metaFrame.channel(fInds(1):fInds(end)),opt.stimterm));
frapInd=frapInd+fInds(1)-1;
ts.frapInSec=metaFrame.timeAfter(frapInd);

% get tk image series timeAfter values and Id the time value for the
% numBefore frame
tAfter=metaFrame.timeAfter(fInds);
ts.tkInSec=tAfter;
tNumBefore=tAfter(numBefore);   
ts.stimInSec=metaFrame.timeAfter(frapInd);
%subtract the numBefore frame from all frames to get the elapsed time
%stamps
tAfterElapsed=tAfter-tNumBefore;
ts.tkTime=tAfterElapsed;
ts.frapTime=metaFrame.timeAfter(frapInd)-tNumBefore;
    
    
end
ts.frameRate=round(mean(diff(ts.tkTime(1:10))),2);



end
