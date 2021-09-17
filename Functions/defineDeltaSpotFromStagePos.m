function deltaSpot = defineDeltaSpotFromStagePos(rawDataroot,dataSubdir,channelNames, varargin)
% created 2020-07-28 GB: This function is designed to read the stage
% positions text file from the centerStim experiments and determine if the
% the 12px frap spot shift had been incorporated into the experiment. THis function
%is necessary because we empirically measured that the frapSpot. 
%hits the cell at ~12 pixels to the right of where we actually image the spot.
%Recently, we added a [-12 0] pixel shift in the experiment to compensate for
% our observed frap location. The shift is
% documented in the text file when a frap line is recorded. Older
% iterations that did not shift for the frapspot have x,y coordinates only.
% THis function will read the text file, if "FRAP" is discovered the
% deltaSpot is set to [0 0] indicating that the frap coordinates do not
% need to be altered. In contrast, if "FRAP" is not detected, the [12 0]
% delta spost will be set, shifting the frap location to our empirically
% measured location.

%%  process varargin inputs
opt.keyword='centroid';
opt.filename='stage_positions.txt';
opt.numfilenames=[];
opt.numbefore=5;
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%% calc number of filenames
numImFiles=findNumImgNamesInFolder(rawDataroot,dataSubdir,channelNames);
if isempty(opt.numfilenames)
    opt.numfilenames=numImFiles;
end
%% Detect frap shifts

%First find a folder with the right number of images in it. Some folders
%are empty.
useFolderFlag=false;
siteCheck=0;
siteCheckMax=length(dataSubdir);
while ~useFolderFlag && siteCheck<siteCheckMax
    
    siteCheck=siteCheck+1;
    clear im
    folder=[rawDataroot filesep dataSubdir{siteCheck}];
    if  boolRegExp(dataSubdir{siteCheck},opt.keyword)
        %~isempty(find(cellfun(@(x) boolRegExp(dataSubdir{siteCheck}, x),opt.keyword),1))
        
        numBefore = findNumImgBeforeFrapStim(folder,channelNames,'printwarning',false); %calculate the number of frames before stimulus
        if any(size(dir([folder '/*.tif' ]),1)) ==1 && numBefore==opt.numbefore
            filenames1=getFilenames(folder,channelNames{1});
            filenames2=getFilenames(folder,channelNames{2});
            fid=fopen([folder filesep opt.filename],'r');
            lines=readAllLines(fid);
            fclose(fid);
            if ~isempty(filenames1) && length(filenames1) == opt.numfilenames && length(lines)>1
                useFolderFlag=true;
            end
        end
    end
    
end
%read text file and determine if there is a frap shift
frapDetected=any(cellfun(@(x) boolRegExp(x,'Frap'),lines));
if frapDetected
    deltaSpot=[0 0];
else
    deltaSpot=[12 0];
end

end
