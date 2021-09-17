function [p, imgF, imAlign] = collectImgAndAlignWithPQuadFit(root, dataSubdir, channelName, varargin)
%% this function will read in image files from a specified data folder,
%and then use Sean's computeRegistrationParametersQuadFit function to
%generate the alignment parameters. The channelName variable should be a
%1x2 cell array that contains both channel names as strings. The optional
%P0 provides the computeRegistrationParametersQuadFit function with an
%initial starting guess for the translation part of the alignment. P0 will
%be useful for TomKat image alignment. 


%2021-05-29: added imAlign as function output. imAlign contains the
%prealginment images. Use this for modifying p0 if the first execution of
%htis function fails.

%2021-05-29: Courtney and Kacey's High mag exp cells have lots of
%background from the RFP cells that confuse alginment algorithm. I added the 
%optional flatfield correction that we use for detecting cells on the
%microscope on the fly. The varargin input is a TF question for using the
%correction.
%% read in the varargins
opt.p0=[];
opt.binsize=1;
opt.filenumlimit=10;
opt.ffctf=false;

imgName1=channelName{1};
imgName2=channelName{2};

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end



 %% Compute an optimal image alignment
 
 %allocate space for  
        numPixels=1024;
        imAlign=zeros(numPixels,numPixels,2);
        temp1=zeros(numPixels,numPixels,length(dataSubdir));
        temp2=zeros(numPixels,numPixels,length(dataSubdir));
        
        count=0; cellCounter=0;
        while count < length(dataSubdir) && cellCounter < opt.filenumlimit
            count=count+1;
            folder=[root filesep dataSubdir{count}];
            filenames1=getFilenames(folder, imgName1);
            filenames2=getFilenames(folder, imgName2);
            if ~isempty(filenames1)
                cellCounter=cellCounter+1;
            temp1(:,:,count)=imageBin(double(imread([folder filesep filenames1{1}])), opt.binsize);
            temp2(:,:,count)=imageBin(double(imread([folder filesep filenames2{1}])), opt.binsize);
            else
                continue
            end
        end
           
        imAlign(:,:,1)=max(temp1,[],3);
        imAlign(:,:,2)=max(temp2,[],3);
        % optional flat field correction if bg noise is interfering with
        % alignment
        if opt.ffctf
            s1=1024;
            s2=1024;
            pFit=[565 522 320];  % using the center from the fitting above, and a manually chosen gaussian width
            [xGrid,yGrid]=meshgrid(1:s2,1:s1);
            imFFC=exp(-1*((xGrid-pFit(1)).^2 + (yGrid-pFit(2)).^2)/(2*pFit(3)*pFit(3)));
            scanImage=(double(imAlign-100)./imFFC);
            scanImage=double(imageSubtractBackground(scanImage,40,256));
            imAlign=scanImage;
        end
        
        p=computeRegistrationParametersQuadFit(imAlign,opt.p0);
        
        imgF=registerImagesFromQuadFit(imAlign,[p.pX p.pY]);
%        imgF=registerImagesFromQuadFit(imAlign,[opt.p0.pX opt.p0.pY]);
        imgF=imageBin(imgF,1);
%         showImagesMergeChannels(imgF(:,:,1),imgF(:,:,2));
        [getErrorValuesForRegisteredImages([p.pX p.pY],imgF) getErrorValuesForRegisteredImages([p.pX p.pY],imgF)] %imstack for getError... was img0 and imgA0. Not sure what these were, I changed to imgF-GB  
end