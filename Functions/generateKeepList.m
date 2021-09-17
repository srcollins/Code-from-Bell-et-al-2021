function [keepList,keepListFoldNames]=generateKeepList(rawDataroot,dataSubdir,channelNames,p,dirSearchTerm,varargin) 
%% This function is a screening tool to allow the user to watch a movie of
%a cell to determine if the cell should be kept for further processing. The
%fuction was upgraded on 1/21/2020 to show the user a movie of the cell
%rather than a single frame. Additionally, the function will check if a
%frap stimulus was applied and scatter a frap spot on the image. 

% 2020-04-29 update to function to add opt.keyword input. This allows the
% user to specify 'keywords' that will allow the function to descriminate
% different dataSubdir folders. My goal is to weed out carole parents cell
% lines that we are not using for my current paper with out changing the
% indicies that we use for identifying the keepList. I also added
% dirSearchTerm as a required function input incase the keyword is left
% blank. The dirSearchTerm is the string used to identify the folders in
% dataSubdir.

% 2020-07-31 update: altered the function to output the foldername in
% addtion to the keeplist position.
opt.frapcoors=[];
opt.numfilenames=[];
opt.cropsize=[];
opt.keyword={};

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% calc number of filenames
numImFiles=findNumImgNamesInFolder(rawDataroot,dataSubdir,channelNames);
if isempty(opt.numfilenames)
    opt.numfilenames=numImFiles;
end
fprintf('NumFolders:%i',length(dataSubdir));
%% if opt.keyword is empty scan datasubdir for common search term
if isempty(opt.keyword)
    opt.keyword={dirSearchTerm};
end
%% define keep list
keepList=[];
keepListFoldNames={};
prompt='keep this folder?';
dataSubdir=natsortfiles(dataSubdir);
for siteCheck=1:length(dataSubdir)
    
    clear im
    folder=[rawDataroot filesep dataSubdir{siteCheck}];
    if ~isempty(find(cellfun(@(x) boolRegExp(dataSubdir{siteCheck}, x),opt.keyword),1))
        
        %numBefore = findNumImgBeforeFrapStim(folder,channelNames);
        %%calculate the number of frames before stimulus SH 8-18-21
        numBefore = 10;
        if any(size(dir([folder '/*.tif' ]),1)) ==1 && ~isempty(numBefore)
            filenames1=getFilenames(folder,channelNames{1});
            filenames2=getFilenames(folder,channelNames{2});
            
            if ~isempty(filenames1) && length(filenames1) == opt.numfilenames
                figure(1);
                for frameNum=1:length(filenames1)
                    im(:,:,1)=double(imread([folder filesep filenames1{frameNum}]));
                    im(:,:,2)=double(imread([folder filesep filenames2{frameNum}]));
                    im2=imageBin(registerImagesFromQuadFit(im,p),1);
                    im2=cropImMidOut(im2,'cropsize', opt.cropsize);
                    showImagesMergeChannels(im2(:,:,1),im2(:,:,2));
                    drawnow; hold on;
                    text(20,20,sprintf('%s',dataSubdir{siteCheck}),'Color','w','FontSize',15);
                    if ~isempty(opt.frapcoors) && frameNum> numBefore
                        hold on; scatter(opt.frapcoors(1), opt.frapcoors(2), 'r', 'LineWidth', 3);
                        %                 scatter(opt.frapcoors(2), opt.frapcoors(1), 'w', 'LineWidth', 3);
                    end
                    pause(0.05);
                end %img play movie loop
                pause(0.5);
                s=input(prompt,'s');
                if boolRegExp(s,'^[YyKk]')
                    keepList=[keepList siteCheck];
                    keepListFoldNames=[keepListFoldNames folder];
                    disp('Cell kept')
                end
                close(1);
                fprintf('wells remaining:%i\n',length(dataSubdir)-siteCheck);
            end %length(filenames1) > 1
        end % if any(size(dir([folder '/*.tif' ]),1)) ==1
    end %if ~isempty(find(cellfun(@(x) boolRegExp(dataSubdir{siteCheck}, x),opt.keyword),1))
end %sitecheck

end
