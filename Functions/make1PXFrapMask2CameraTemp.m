function [frapMask1PX, frapInd, frapMask1PXrt,frapIndrt] = make1PXFrapMask2CameraTemp(root, p, varargin)
% Function inputs: root, p-paramete,folder 'keyword', imageID 'keyWord', crop Range of final Im (ie [600 600]), translate
% spot ([x,y]; Use name value pairs with names coming from the opt
% structure
% function returns a logical mask of a 1 pixel frapspot based on the max
% intensity measurement of the actual frap image. 

% The function now returns the frapInd in [X Y] order as of 1/8/2020; The
% older versions of this function returned the frap spot in row by column
% format which was [y,x]; I switched it because I am just using the frapInd
% for scattering the frap spot onto images and the scatter function requies
% the [X Y ] format.


%The translate spot parameter will shift the spot by the specified pixel distances
% eg [12, 0] = 12 px to the right
% NOTE this function is designed for fret pair experiments and requires a p
% parameter
%% varargin inputs
tempRoot=root;

opt.deltaspot=[0 0];
opt.cropsize=[];
opt.imid='FRAP';
opt.keyword='FRAP';



for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

if nargin < 2
    fprintf('data root and alignment parameters are required.');
    return
end


%% Find the frap Image in the folder
    clear frapTemp;
    
% first Check if frap Im is in the root provided
fileNameCheck=getFilenames(tempRoot,opt.imid);
if ~isempty(fileNameCheck)
  filenames1=getFilenames(tempRoot,opt.imid);
  folder=tempRoot;
else
    % set root dir and find sub folders with FRAP as key word
    frapDir=getSubdirectories(tempRoot,opt.keyword);
    if isempty(frapDir)
        idcs=strfind(tempRoot,filesep);
        tempRoot = tempRoot(1:idcs(end)-1);
        frapDir=getSubdirectories(tempRoot,opt.keyword);
        frapDirChar=char(frapDir{1});
    else
        frapDirChar=char(frapDir{1});
    end
    if isempty(frapDirChar)
    sprintf('please define a case sensitive keyword to find the correct sub directory');
        return
    end
    %define folder and filenames
    folder=[tempRoot filesep frapDirChar];
    filenames1=getFilenames(folder,opt.imid);
    if isempty(filenames1)
        sprintf('please define a case sensitive keyword to find the correct filenames');
        return
    end    
end %Check file Names in TempRoot.
    %%
    frapIm=[];
    for i=1:length(filenames1)
        frapIm(:,:,i)=double(imread([folder filesep filenames1{i}]));
    end
    if size(frapIm,3)>1
        frapIm=sum(frapIm,3);
    end
    
    frapTemp(:,:,1)=(frapIm);
    frapTemp(:,:,2)=(frapIm);
    frapTemp=registerImagesFromQuadFit(frapTemp,p);    
    frapTemp=cropImMidOut(frapTemp,'cropsize', opt.cropsize);
    
    [~,frapCenter]=max(vect(frapTemp(:,:,1)));
    
    [fInd(1),fInd(2)]=ind2sub(size(frapTemp(:,:,1)),frapCenter); %(ind2sub is a row by column format.
    % make the non-shifted mask
    tempMask=false(size(frapTemp(:,:,1)));
    tempMask(fInd(1),fInd(2))=1; %positioning the mask is RbyC so 
    frapMask1PX=tempMask;
    frapInd=convertXYtoRC(fInd);
    
    % incorporate the deltaSpot into frapInd. fInd is stll in Row by Col
    % notation while the delta spot is XY convention.
    fInd(2)=fInd(2)+opt.deltaspot(1);
    fInd(1)=fInd(1)+opt.deltaspot(2);
    tempMask=false(size(frapTemp(:,:,1)));
    tempMask(fInd(1),fInd(2))=1;
    frapMask1PXrt=tempMask;
    frapIndrt=convertXYtoRC(fInd);
   
end
    
    