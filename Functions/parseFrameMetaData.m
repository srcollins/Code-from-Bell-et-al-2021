function metaFrame = parseFrameMetaData(dataDir, varargin)
% This function is desinged to read all metadata files from an input
% folder. THe funciton will combine the data structures and eliminate empty
% nan values. This ouput info can then be used to find the timing for each
% cell frame capture


%2020-07-30 Update: found a metaframe that had duplicate images with
%duplicate time stamps. I am not sure how this happened but I will remove
%duplicates in this function. *) I found the error. A metaData file was
%created but not used, thus it lacked frame and event variables, throwing
%off my function. Now the function checks to make sure frame exists before
%proceeding.

%2021-02-24 UV uncaging is getting counted in metaFrame causing duplicates.
% I will have the function find the uncs and remove them.

%% process varargin

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% 
metaFiles=getFilenames(dataDir,'metadata');
if isempty(metaFiles)
    metaFrame=[];
    return
end

% load the first metadata file to get fieldNames
load([dataDir filesep metaFiles{1}]);
metaFrame=frame;
fNames=fieldnames(metaFrame);
fRange=[2:length(metaFiles)];

% combine all metaFiles into metaFrame a single struct array
for i=1:length(metaFiles)-1
    clear frame;
    load([dataDir filesep metaFiles{fRange(i)}]);
    if exist('frame') ==1
        for f=1:length(fNames)
            metaFrame.(fNames{f})=[metaFrame.(fNames{f}); frame.(fNames{f})];
        end
    end
end

% remove empty positions.
eInd=isnan(metaFrame.frameID);
for e=1:length(fNames)
    if boolRegExp(fNames{e},'xy')
    metaFrame.(fNames{e})(eInd,:)=[];
    else
         metaFrame.(fNames{e})(eInd)=[];
    end
end

 % remove uncaging lightpulses that end up as duplicates 
uInd=find(boolRegExp(metaFrame.label,'unc'));
for u=1:length(fNames)
    if boolRegExp(fNames{u},'xy')
    metaFrame.(fNames{u})(uInd,:)=[];
    else
         metaFrame.(fNames{u})(uInd)=[];
    end
end
 

end


