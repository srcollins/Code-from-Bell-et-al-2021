function compositeFrapIm = makeCompositeFrapIm(root, p, varargin)
%% This function is designed to build a summed frap image from several frap images 
%in the case where a bright frap image was not collected during the
%experiment. THe summed frap image will be processed seperately using the
%make1PXFrapMask2Camera function. The function can handle a keyword that 
%specifies multiple subdirectories.

%the function will detect if some of the folders are empty and
%automatically skip those

tempRoot=root;

opt.cropsize=[];
opt.imid='FRAP';
opt.keyword='centroid';


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

if nargin < 2
    fprintf('data root and alignment parameters are required.');
    return
end

%% Find the frap Image in the folder
    clear frapTemp;
    % set root dir and find sub folders with FRAP as key word
    frapDir=getSubdirectories(tempRoot,opt.keyword);
    
 
    if isempty(frapDir)
    sprintf('please define a case sensitive keyword to find the correct sub directory');
        return
    end
   
  %% read in frap images and return a summed frap image  
  %define folder and filenames
  xCent=512;
  yCent=512;
  fCount=0;
  goodFrapCount=0;
  for foldNum=1:length(frapDir)
      clear frapTemp
      folder=[tempRoot filesep frapDir{foldNum}];
      filenames1=getFilenames(folder,opt.imid);
      
      %only take folders with frap filenames
      if ~isempty(filenames1)
          fCount=fCount+1;
          goodFTemp=nan(1024,1024,length(filenames1));
          
          for i=1:length(filenames1) %read in frap Images
              tempFrapIm=double(imread([folder filesep filenames1{i}]));
              imagesc(tempFrapIm,[100 300]);
              
              %detect max pixle intensity and only keep image if max is
              %near the center
              [~,ind]=max(vect(tempFrapIm));
              [x,y]=ind2sub([1024,1024],ind);
              hold on; scatter(x,y,'r','LineWidth',2);
              pause(0.3);
              xCheck=x>(xCent-100) & x<(xCent+100);
              yCheck=y>(yCent-100) & y<(yCent+100);
              if xCheck && yCheck
                  goodFrapCount=goodFrapCount+1;
                  goodFTemp(:,:,goodFrapCount)=tempFrapIm;
              end
              
          end
      else
          continue %skip this iteration of the loop if there are no frap images
      end
     
      frapIm=nansum(goodFTemp,3); %sum the frap Images to make one image
      frapTemp(:,:,1)=(frapIm);
      frapTemp(:,:,2)=(frapIm);
      frapTemp=registerImagesFromQuadFit(frapTemp,p); %register the image
      frapTemp=cropImMidOut(frapTemp,'cropsize', opt.cropsize); %crop image
      frapTempAll(:,:,fCount)=frapTemp(:,:,1); %store the frap im for each data subDir;
  end
  
%    if isempty(filenames1)
%           sprintf('please define a case sensitive keyword to find the correct filenames');
%           return
%       end
  compositeFrapIm=nansum(frapTempAll,3);