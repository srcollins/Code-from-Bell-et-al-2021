function inds2Plot = sortpDatPLBbyStimConds(pDatPDL,cellType,Power,numStim,varargin)
%2020-08-07: This function will take key words pertaining to celltype,
%stimpower and the number of stimulations and search the stimGroups stored
%in pDatPDL. THe function will then return the specific indicies of the
%pDatPDL array to plot. Input 'all' for any of the keyword variables if you
%want the function to return all groups for that key word.

%eg cellType='all'; power=5000; numstim=30; will return all cells with the
%matching power and numstim vals

% make the input for any of the conds a cell array for multiple conditions
%% process varargin inputs


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% Find the stimInds

cellTinds=[];
pwrInds=[];
numStimInds=[];
% sort keywords: cellType
if any(boolRegExp(cellType,'all'))
    cellTinds=find(true(1,length(pDatPDL)));
else
    for c=1:length(cellType)
    tempCinds=find(arrayfun(@(x) boolRegExp(x.tCond,cellType{c}),pDatPDL));
    cellTinds=[cellTinds tempCinds];
    end
end
%sort keywords: stime Power
if any(boolRegExp(Power,'all'))
    pwrInds=find(true(1,length(pDatPDL)));
else
    for p=1:length(Power)
    tempPinds=find(arrayfun(@(x) boolRegExp(x.pwr,Power{p}),pDatPDL));
    pwrInds=[pwrInds tempPinds];
    end
end

%sort keywords: number of stims
if any(boolRegExp(numStim,'all'))
    numStimInds=find(true(1,length(pDatPDL)));
else
    for n=1:length(numStim)
    tempNinds=find(arrayfun(@(x) boolRegExp(x.numStim,numStim{n}),pDatPDL));
    numStimInds=[numStimInds tempNinds];
    end
end

inds2Plot=intersect3(cellTinds,pwrInds,numStimInds);
end