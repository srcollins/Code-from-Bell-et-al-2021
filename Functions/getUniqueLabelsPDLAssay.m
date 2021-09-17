function cLabel=getUniqueLabelsPDLAssay(dataSubdir,cellLine,varargin)
% 10/09/20 GB: This function is designed to find cell condition
% labels for the PDL assay. cLabel will be a data structure that will
% contain the cell treatment label, the number of stimulations and the
% light stimulus power.

%% process varargins
opt.wellpatt='(^[A-Z][0-9][0-9])';
opt.int='(?<=Int_)[0-9]{0,3}';
opt.numstim='(?<=numStims_)[0-9]{0,3}';
opt.exp='(?<=Exp_)[0-9]{0,3}';
opt.constructlabels='SC[0-9]{0,4}';
opt.keepplasmidlabels=false;
opt.sharedconditions=''; %input characters to be removed from cLabel.tCond to make tCondUnique (ex. plasmids or treatment shared by all cells)
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% pre allocate label

for i=length(dataSubdir):-1:1
    cLabel(i).tCond=[];
    cLabel(i).tCondUnique=[];
    cLabel(i).numStim=[];
    cLabel(i).pwr=[];
    cLabel(i).dataSubdir=[];
    cLabel(i).stimGroup=[];
    
end
%% extract dataSubDir info and build cLabel
cellCondPatt=sprintf('(?<=%s_|-).*(?=_Int)',cellLine);

for d=1:length(dataSubdir)
    cLabel(d).dataSubdir=dataSubdir{d};
    % define treatment condition for the cells
    cellCond=regexpi(dataSubdir{d},cellCondPatt,'match'); % get contents of fold label between 'PLB' and '_Int'
    % id plasmid constructs and remove these labels from cellCond
    dnaTemp=regexp(cellCond,opt.constructlabels,'match');
    plasmidPat=sprintf('(%s.*%s)',dnaTemp{1}{1},dnaTemp{1}{end});
    plasmidLabel=regexp(cellCond,plasmidPat,'match');
    if opt.keepplasmidlabels == false
        cellCond=regexprep(cellCond,plasmidLabel{1},''); % remove plasmid constructs from cellCond
    end
    if boolRegExp(cellCond,'^(-|_)')
        cellCond{1}(1)=[];
    end
    if boolRegExp(cellCond,'(-|_)$')
        cellCond{1}(end)=[];
    end
    cLabel(d).tCond=cellCond{1};
    
    %remove tCond shared by all cells to make tCondUnique -default makes an
    %identical label. String removed specified by opt.sharedconditions. 
    cellCondUnique = regexprep(cellCond,opt.sharedconditions,'');
    cellCondUnique = regexprep(cellCondUnique,'--','-');
    cLabel(d).tCondUnique = cellCondUnique{1};
    
    %define the number of stims
    numStim=regexp(dataSubdir{d},opt.numstim,'match');
    if isempty(numStim)
        numStim={'1'};
    end
    cLabel(d).numStim=numStim{1};
    % define stim 1ight power: Int*Exp=Pwr
    int=regexp(dataSubdir{d},opt.int,'match');
    exp=regexp(dataSubdir{d},opt.exp,'match');
    pwr=str2double(int)*str2double(exp);
    cLabel(d).pwr=num2str(pwr);
    cLabel(d).stimGroup=sprintf('%s-%s-%s',cLabel(d).tCond,cLabel(d).pwr,cLabel(d).numStim);
end
