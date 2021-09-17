function deleteInds=removeCellsUsingmovieKeepLists(dir,bigRes,masterRoot,varargin)

%% 2020-06-05 GB: Remove specific cell data from bigRes.
%After watching the movies that got created along with this data, I found
%some additional cells that need to be removed. I will use the fileName to
%Id the potion in the structure and then clean out the data for that cell.
%This function just Identifies the indicies for the cells that need to be
%removed

%% process Varargin

%% process varargin inputs
opt.fnamekeyword='Ret';


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% find indicies for cells to remove from bigRes.

%Build array counter
    stimLabelsAll=arrayfun(@(x) x.name{1},bigRes,'UniformOutput',false);
    deleteInds=cell(1,length(stimLabelsAll));
    
    
    for d=1:length(dir)
        parentFolder=dir{d};
        scriptR=[masterRoot 'Scripts' filesep parentFolder];
        saveR=[masterRoot 'Analyzed Data' filesep parentFolder];
        load([scriptR filesep 'movieKeepLists.mat']);
        
        
        % I found that there are duplicate fileNames in bigRes from different
        % experiments. I also need to use the parentFoldre name to make sure I
        % am deleting the right cell.
        % newer files have Ret in teh file name. Dont look at the older ones
    fileNameInd=boolRegExp(notKeepersFileNames,opt.fnamekeyword);
    notKeepersFileNames=notKeepersFileNames(fileNameInd);
        
        
        for n=1:length(notKeepersFileNames)
            tempStr=notKeepersFileNames{n};
            tempLabelType=regexp(tempStr,'(?<=PLB[-_]).*(?=Cell)','match');
            bigRInd=find(strcmpi(stimLabelsAll,tempLabelType));
            badCellInd=find(boolRegExp(bigRes(bigRInd).fileName,[tempStr '$']));
            
            if ~isempty(badCellInd)
                for bc=1:length(badCellInd)
                    if boolRegExp(bigRes(bigRInd).parentFolder{badCellInd(bc)},parentFolder)
                        deleteInds{bigRInd}=[deleteInds{bigRInd} badCellInd(bc)];
                        
                    end %if ~isempty(badCellInd)
                end % for bc
            end %for n=1:length(notKeepersFileNames)
        end %for d
    end