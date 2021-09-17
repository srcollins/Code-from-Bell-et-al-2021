function failFlag=findCellsWithLgSqMask(fNamePath,varargin)
% when creating imDat for the centerStim Experiments, I would make the mask
% a huge square if the mask image was empty. This was done to prevent
% errors, however, the strategy wasnt't great and is now causing errors
% else where. Specifically, these square masks are making it through to all
% of the data analysis steps including the diffExp plots and the
% kymographs.

% This function will read in all imDat files based on the indicies for
% res. It will then check if any of the frames in the range that we care
% about (2:16) have this square mask. If they do, the function will save
% save fail flag as true.


%% process varargin inputs

opt.framerange=[2:20];


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%% read in imDats and find bad cells
load(fNamePath,'imDat');
t=imDat.cropMask;
fRange=[1:opt.framerange(end)];
t=t(fRange);  
sqTestMask=zeros(size(t{1}));
sqTestMask(5:size(sqTestMask,1)-5, 5:size(sqTestMask,2)-5)=1;

maskComp=cellfun(@(x) all(x==sqTestMask,'all'), t,'UniformOutput',false);
maskCompVect=cellfun(@(x) x,maskComp);

failFlag=any(maskCompVect(opt.framerange));
end
