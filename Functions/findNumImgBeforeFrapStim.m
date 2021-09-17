function numBefore = findNumImgBeforeFrapStim(folder,imgChannelNames, varargin)
% This function will read all of the filenames in an image folder and count
% the number of frames before the frap stimulus. In many cases, we are
% imaging the tirf channels before the frap stimulus, thus if we set the
% numBefore to 5 then the stimulus would actually be applied after the
% first fret pair were imaged for the num after. To compensate, I added the
% optional imbeforestim variable. If imgbeforestim=true,
% numbefore=numbefore+1;
opt.imgbeforestim=false;
opt.printwarning=true;
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%%
% load fileNames
allImg=getFilenames(folder,'img');
frapInd=find(boolRegExp(allImg,'FRAP'),1);
if isempty(frapInd) && opt.printwarning
    fprintf('Frap Image was not discovered, numBefore set to empty vector.\n');
    numBefore=[];
    return
end
% I am updating how I count the number of frames before as my current
% strategy is prone to errors
tirfInd=find(boolRegExp(allImg,imgChannelNames{1}));
tirfLogicals=tirfInd<frapInd;
numBefore=sum(tirfLogicals);

%old method, replaced 2020/02/14
% tirfInd=find(boolRegExp(allImg,imgChannelNames{1}));
% doubleNumImg=max(tirfInd(tirfInd<frapInd));
% numBefore=doubleNumImg/2;
if opt.imgbeforestim 
    numBefore=numBefore+1;
end

end