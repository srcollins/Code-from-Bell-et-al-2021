function compositeBleachCorr=buildCompositeBleachCorrectionForMultiStim(bleachCorrectionImage,imDat,varargin)
%2021-02-18 GB: This function builds a composite photobleach correction 
%image for the centerstim experiments where multiple sequential
%frap stimulations were delivered. I need to double check with Sean,
%but here is my plan. I will give the function single bleach correction that has been
%calculated by makeLocalBleachCorrectionImage and the imDat.times struct, which contains
%the number of frap stimulations. The correction essentially needs to be 
%refreshed with every new stimulation, eg the correction for frapStim2
%(frame12)=bleachCorrImg{12}.*bleachCorrImg{11}
%% process varargin inputs

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% convert bleachCorrection Image to 3d matrix
bleachMat=nan(size(bleachCorrectionImage{1},1),size(bleachCorrectionImage{1},2),length(bleachCorrectionImage));
compositeBleachCorr=cell(1,length(bleachCorrectionImage));
for L=1:length(bleachCorrectionImage)
    bleachMat(:,:,L)=bleachCorrectionImage{L};
end

%% detect number of frap stimulations and build bleachImg mat for each frap iteration.
lenFrap=length(imDat.times.frapTime);
for i=lenFrap:-1:1
    b(i).corr=[];
end
onesImg=bleachCorrectionImage{1};
b(1).corr=bleachMat;
for i=2:lenFrap
 b(i).corr=b(i-1).corr;
 b(i).corr(:,:,end)=[];
 b(i).corr=cat(3,onesImg, b(i).corr);   
end

%% build final composite correction
for c=1:length(compositeBleachCorr)
    tempIm=nan(size(bleachCorrectionImage{1},1),size(bleachCorrectionImage{1},2),length(b));
    for i=1:length(b)
    tempIm(:,:,i)=b(i).corr(:,:,c);
    end
   compositeBleachCorr{c}=prod(tempIm,3);
end
%compositeBleachCorr
