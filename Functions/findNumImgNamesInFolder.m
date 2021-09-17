function numImFiles=findNumImgNamesInFolder(rawDataroot,dataSubdir,channelNames)
%This function will look through a number of exp folders, record the number
%of file names for each img
numF1=nan(1,length(dataSubdir));
numF2=nan(1,length(dataSubdir));

for i=1:length(dataSubdir)
    folder=[rawDataroot filesep dataSubdir{i}];
    filenames1=getFilenames(folder,channelNames{1});
    filenames2=getFilenames(folder,channelNames{2});
    numF1(i)=length(filenames1);
    numF2(i)=length(filenames2);
end
% if too many zeros, that will kill the mode

%filter out zeros
numF1=numF1(numF1>0);
numF2=numF2(numF2>0);
modeF1=mode(numF1);
modeF2=mode(numF2);
if modeF1==modeF2
    numImFiles=modeF1;
else
    fprintf('length Filenames is not equal\n');
    return
end

end
