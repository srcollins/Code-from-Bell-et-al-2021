function cLabelOut=cLabelAdjPDLassay(cLabel,searchStr,replaceStr,varargin)
% 2020-10-17 GB: THis function is designed to read the data sturcture
% output from getUniqueLabelsPDLAssay and adjust the tCond label names to
% remove unique but similar labels. 

%inputs:
%1) cLabel data structure
%2) searchStr cell array: contains strings that will be used to id specific
%labels that need adjusting (eg 'Ctrl', 'Ctrl-Ret')
%3) replaceStr cell array: must be same length as searchStr. This cell
%array contains the labels that will be used to replace the tConds
%identified by the searchStr;

%% process varargin
opt.all=false;
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%% check lengths of searchStr and replaceStr

if numel(searchStr)~=numel(replaceStr)
    fprintf('searchStr and replaceStr must be equal length');
    return
end
%% define cLabelOut
cLabelOut=cLabel;


%% update cLabelout.tCond and cLabelOut.stimGroup
if opt.all
    for i=1:length(cLabelOut)
        cLabelOut(i).tCond=replaceStr{1};
        cLabelOut(i).stimGroup=sprintf('%s-%s-%s',cLabelOut(i).tCond,cLabelOut(i).pwr,cLabelOut(i).numStim);
    end
    
else
    
    for s=1:length(searchStr)
        ind=[];
        ind=find(arrayfun(@(x) boolRegExp(x.tCond,searchStr{s}),cLabelOut));
        
        for i=ind
            cLabelOut(i).tCond=replaceStr{s};
            cLabelOut(i).stimGroup=sprintf('%s-%s-%s',cLabelOut(i).tCond,cLabelOut(i).pwr,cLabelOut(i).numStim);
        end
        
    end
end