function resOut=cleanUpBigRes(bigRes,masterRoot, varargin)
%% 2021-02-19 GB: THis code got to be too much to look at. Its designed to remove
%duplicate structures and then duplicate cells from res struct arrays that
%were concatinated with duplicates.

%% process varargin inputs

opt.deleteduplicatecelldataTF=false;
opt.getdirkeyword='PP1';

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end


%% clean up bigRes
%When analyzing the data I had to interupt the analysis part way
%through...resulting in two structures that need to be concatonated. BigRes
% has a few empty values

emptyInd=find(arrayfun(@(x) isempty(x.name),bigRes));
bigRes(emptyInd)=[];
% find overlapping 
uniqueLabels=unique(arrayfun(@(x) x.name(1), bigRes));

for b=1:length(bigRes)
    labelInd(b)=find(strcmp(uniqueLabels,bigRes(b).name{1}));
end

% Here is the painful part, duplicates cant simply be concatenated unless
% they are the same size, so I need to find the dups and then concatenate
% on a per field basis. 
dupCounter=0;
for i=1:length(bigRes)
    
    dups=find(labelInd==labelInd(i));
    if length(dups)>1
        
        for j=2:length(dups)
            dupCounter=dupCounter+1;
        dup2delete(dupCounter)=dups(j);
        %updated 2020-08-04 for all new metrics
        bigRes(dups(1)).name=[bigRes(dups(1)).name bigRes(dups(j)).name];
        bigRes(dups(1)).dist=[bigRes(dups(1)).dist bigRes(dups(j)).dist];
        
      
        %ref5
        bigRes(dups(1)).ref5DiffExp=cellfun(@(x,y) [x ; y], bigRes(dups(1)).ref5DiffExp, bigRes(dups(j)).ref5DiffExp,'UniformOutput', false);
        bigRes(dups(1)).ref5NumExp=cellfun(@(x,y) [x ; y], bigRes(dups(1)).ref5NumExp, bigRes(dups(j)).ref5NumExp,'UniformOutput', false);
        %ref2-5
        bigRes(dups(1)).ref2_5DiffExp=cellfun(@(x,y) [x ; y], bigRes(dups(1)).ref2_5DiffExp, bigRes(dups(j)).ref2_5DiffExp,'UniformOutput', false);
        bigRes(dups(1)).ref2_5NumExp=cellfun(@(x,y) [x ; y], bigRes(dups(1)).ref2_5NumExp, bigRes(dups(j)).ref2_5NumExp,'UniformOutput', false);
        %ref2
        bigRes(dups(1)).ref2DiffExp=cellfun(@(x,y) [x ; y], bigRes(dups(1)).ref2DiffExp, bigRes(dups(j)).ref2DiffExp,'UniformOutput', false);
        bigRes(dups(1)).ref2NumExp=cellfun(@(x,y) [x ; y], bigRes(dups(1)).ref2NumExp, bigRes(dups(j)).ref2NumExp,'UniformOutput', false);
        bigRes(dups(1)).parentFolder=[bigRes(dups(1)).parentFolder bigRes(dups(j)).parentFolder];
        bigRes(dups(1)).date=[bigRes(dups(1)).date bigRes(dups(j)).date];
        bigRes(dups(1)).fileName=[bigRes(dups(1)).fileName bigRes(dups(j)).fileName];
        bigRes(dups(1)).numBefore=[bigRes(dups(1)).numBefore bigRes(dups(j)).numBefore];
        % new metrics
        bigRes(dups(1)).wholeCellMeanFret=[bigRes(dups(1)).wholeCellMeanFret bigRes(dups(j)).wholeCellMeanFret];
        bigRes(dups(1)).frap2edge=[bigRes(dups(1)).frap2edge bigRes(dups(j)).frap2edge];
        bigRes(dups(1)).andMaskChange=[bigRes(dups(1)).andMaskChange bigRes(dups(j)).andMaskChange];
        bigRes(dups(1)).euclidDist=[bigRes(dups(1)).euclidDist bigRes(dups(j)).euclidDist];
        bigRes(dups(1)).cellCentroids=[bigRes(dups(1)).cellCentroids bigRes(dups(j)).cellCentroids];
        
        
        end %j
    end
end

dup2delete=unique(dup2delete);
bigRes(dup2delete)=[];

%% actually delete the bad cells from bigRes
if opt.deleteduplicatecelldataTF
    dir=getDirectories([masterRoot 'Analyzed Data'],opt.getdirkeyword );
    deleteInds=removeCellsUsingmovieKeepLists(dir,bigRes,masterRoot);
    for r=1:length(bigRes)
        
        bigRes(r).name(deleteInds{r})=[];
        bigRes(r).date(deleteInds{r})=[];
        
        bigRes(r).parentFolder(deleteInds{r})=[];
        bigRes(r).fileName(deleteInds{r})=[];
        bigRes(r).numBefore(deleteInds{r})=[];
        
        bigRes(r).frap2edge(deleteInds{r})=[];
        bigRes(r).andMaskChange(deleteInds{r})=[];
        bigRes(r).euclidDist(deleteInds{r})=[];
        bigRes(r).cellCentroids(deleteInds{r})=[];
        
        bigRes(r).wholeCellMeanFret(deleteInds{r})=[];
        
        
        for dif=1:length(bigRes(r).ref5DiffExp)
            bigRes(r).ref5DiffExp{dif}(deleteInds{r},:)=[];
            bigRes(r).ref5NumExp{dif}(deleteInds{r},:)=[];
            bigRes(r).ref2_5DiffExp{dif}(deleteInds{r},:)=[];
            bigRes(r).ref2_5NumExp{dif}(deleteInds{r},:)=[];
            bigRes(r).ref2DiffExp{dif}(deleteInds{r},:)=[];
            bigRes(r).ref2NumExp{dif}(deleteInds{r},:)=[];
            
        end %for dif
        
    end
end
%% find duplicate cellData. 
% It looks like some cells may have up to tripicate copies. I need to hunt
% out the triplicates and remove them.

for b=1:length(bigRes)
    [C, ia,~]=unique(bigRes(b).fileName);
    dupInd=1:length(bigRes(b).fileName);
    dupInd(ia)=[];
    
    if length(C)~=length(bigRes(b).fileName)
     
        for c=1:length(C)
            specificDupInd=find(boolRegExp(bigRes(b).fileName,[C{c} '$']));
            uniqueFolderNames=unique(bigRes(b).parentFolder(specificDupInd));
            pFoldInd=cellfun(@(x) find(boolRegExp(bigRes(b).parentFolder(specificDupInd),x)), uniqueFolderNames, 'UniformOutput', false);
            deleteIndsAll=[];
            for pf=1:length(pFoldInd)
              deleteInds=specificDupInd(pFoldInd{pf}(2:end));  
                    deleteIndsAll=[deleteIndsAll deleteInds];
            end% for pf 
                   bigRes(b).name(deleteIndsAll)=[];
                   bigRes(b).date(deleteIndsAll)=[];
                   bigRes(b).parentFolder(deleteIndsAll)=[];
                   bigRes(b).fileName(deleteIndsAll)=[];
                   bigRes(b).numBefore(deleteIndsAll)=[];
                   
                   bigRes(b).wholeCellMeanFret(deleteIndsAll)=[];
                   bigRes(b).frap2edge(deleteIndsAll)=[];
                   bigRes(b).andMaskChange(deleteIndsAll)=[];
                   bigRes(b).euclidDist(deleteIndsAll)=[];
                   bigRes(b).cellCentroids(deleteIndsAll)=[];
                   
                   for dif=1:length(bigRes(b).ref5DiffExp)
                       
                       
                       bigRes(b).ref5DiffExp{dif}(deleteIndsAll,:)=[];
                       bigRes(b).ref5NumExp{dif}(deleteIndsAll,:)=[];
                       bigRes(b).ref2_5DiffExp{dif}(deleteIndsAll,:)=[];
                       bigRes(b).ref2_5NumExp{dif}(deleteIndsAll,:)=[];
                       bigRes(b).ref2DiffExp{dif}(deleteIndsAll,:)=[];
                       bigRes(b).ref2NumExp{dif}(deleteIndsAll,:)=[];
                       
                   end %for dif
            
        end %for c 
    end %if length(C)~=length(bigRes(b).fileName)
end % for length bigRes

%% filter data based on change in AND mask area, mask centroid distances and chage in frap spot position to cell edge
% Sean suggested that I just remove the cell from all analysis rather than
% drop specific time points for a cell. Its more justifiable this way
% unless we are removing lots of cells.

removalCount=zeros(1,length(bigRes));
removalInd=cell(1,length(bigRes)); 

for b=1:length(bigRes)
    for c=1:length(bigRes(b).fileName)
        imDName=sprintf('%s-imgData.mat',bigRes(b).fileName{c});
        fNamePath=[masterRoot 'Analyzed Data' filesep bigRes(b).parentFolder{c} filesep imDName];
        failFlag=findCellsWithLgSqMask(fNamePath,'framerange',[2:20]);
        f2eInd=find(bigRes(b).frap2edge{c}<=4);
        andMInd=find(bigRes(b).andMaskChange{c}>=35);
        euclidDInd=find(bigRes(b).euclidDist{c}>=12)'; %~4um
        if ~isempty(f2eInd) && ~isempty(andMInd) && ~isempty(euclidDInd) || failFlag
            intersectInd=intersect3(f2eInd, andMInd,euclidDInd);
           
            if ~isempty(intersectInd) || failFlag
                removalCount(b)=removalCount(b)+1;
                removalInd{b}=[removalInd{b} c];
                
%                 for i=1:length(intersectInd)
%                     bigRes(b).diffExp{intersectInd(i)}(c,:)=nan;
%                     bigRes(b).numExp{intersectInd(i)}(c,:)=nan;
%                 end
            end
        end
    end
end

%% make a table to show Sean
tNames=arrayfun(@(x) x.name{1}, bigRes,'UniformOutput', false);
treatNames=tNames;

% treatNames=cell(1,length(tNames));
% for i=1:length(tNames)
%     treatNames{i}=tNames{i}{1};
% end

removalCount=removalCount';

numCellsPreRemoval=arrayfun(@(x) length(x.fileName),bigRes)';
cellsRemaining=numCellsPreRemoval-removalCount;

table(removalCount,numCellsPreRemoval,cellsRemaining,'RowNames',treatNames)
%% Delete the cells that are excluded by the 3 metrics in the block above
for b=1:length(bigRes)
    
    bigRes(b).name(removalInd{b})=[];
    bigRes(b).date(removalInd{b})=[];
    bigRes(b).parentFolder(removalInd{b})=[];
    bigRes(b).fileName(removalInd{b})=[];
    bigRes(b).numBefore(removalInd{b})=[];
    
    
    bigRes(b).wholeCellMeanFret(removalInd{b})=[];
    bigRes(b).cellCentroids(removalInd{b})=[];
    bigRes(b).frap2edge(removalInd{b})=[];
    bigRes(b).andMaskChange(removalInd{b})=[];
    bigRes(b).euclidDist(removalInd{b})=[];
    bigRes(b).times(removalInd{b})=[];
    
    for dif=1:length(bigRes(b).ref5DiffExp)
        
        
        bigRes(b).ref5DiffExp{dif}(removalInd{b},:)=[];
        bigRes(b).ref5NumExp{dif}(removalInd{b},:)=[];
        bigRes(b).ref2_5DiffExp{dif}(removalInd{b},:)=[];
        bigRes(b).ref2_5NumExp{dif}(removalInd{b},:)=[];
        bigRes(b).ref2DiffExp{dif}(removalInd{b},:)=[];
        bigRes(b).ref2NumExp{dif}(removalInd{b},:)=[];
        
    end %for dif
    
end


end