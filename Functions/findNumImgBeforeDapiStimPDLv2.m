function numBefore = findNumImgBeforeDapiStimPDLv2(rawDataroot,dataSubdir,imgChannelNames, varargin)
% This function will read all of the filenames in an image folder and count
% the number of frames before the frap stimulus. In many cases, we are
% imaging the tirf channels before the Dapi stimulus, thus if we set the
% numBefore to 5 then the stimulus would actually be applied after the
% first fret pair were imaged for the num after. To compensate, I added the
% optional imbeforestim variable. If imgbeforestim=true,
% numbefore=numbefore+1;

%2020-10-15 Update: SOme of the folders will not have a stimulus applied as
%these are control files. Instead, I will alter the function to read all of
%the folders in the dataSubdir, and find the numbefore position for the
%folders that do have a stimulus. This function will need to be moved
%outside of the dataSubdir for loop.
%% process varargin
opt.stimlabel={'DAPI','UV'};

opt.imgbeforestim=false;
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% preallocate
allNumBefore=nan(1,length(dataSubdir));

%%
for d=1:length(dataSubdir)
    folder=[rawDataroot filesep dataSubdir{d}];
    % load fileNames
    allImg=getFilenames(folder,'img');
    
    %determine stimlabel
    numStims=[0 0];
    for i=1:length(opt.stimlabel)
        numStims(i)=length(find(boolRegExp(allImg,opt.stimlabel{i})));
    end
    
    if any(numStims>0)
        stimLabelInd=numStims>0;
        stimInd=find(boolRegExp(allImg,opt.stimlabel{stimLabelInd}),1);
        tirfInd=find(boolRegExp(allImg,imgChannelNames{1}));
        doubleNumImg=max(tirfInd(tirfInd<stimInd));
        numBefore=ceil(doubleNumImg/2);
        if opt.imgbeforestim
            numBefore=numBefore+1;
        end
    else
        numBefore=0;
    end
    
    allNumBefore(d)=numBefore;
end %d
allNumBefore=allNumBefore(allNumBefore>0);
numBefore=mode(allNumBefore);

end