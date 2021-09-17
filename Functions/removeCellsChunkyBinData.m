function [ratioDat,numRm,numCellTot]=removeCellsChunkyBinData(ratioDat,varargin)
%2020-06-19 GB: I am compiling the chunkybin data into a data structure
%called ratioDat. In ratioDat each chunky bin will contain an nx20 matrix
%where each row represnts a cell and each column is the frame. Cells that
%have at least one nan mean value or a bin area < the area of bin1 are set
%to nan for that bin. This function will search all bins and Id cells that
%failed in one bin and thus should be removed. 

%inputs: 1) ratioDat structure that contains the chunky bin data.
%2) Varargin. Bins closest to the cell edge might have few cell pixels that
%could exclude an otherwise good cell. I will give the user the option to
%determine how many bins to consider. 



%% process varargin inputs
opt.binrange=[1:10];
opt.binname='ratioBin';
opt.notbin0=true;

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% find cells in ratio dat that fail and set them to nan


for i=1:length(ratioDat)
    clear fieldN;
    % determine the number of bins and pull out fieldnames
    fieldNAll=fieldnames(ratioDat(i));
    binInds=find(boolRegExp(fieldNAll,opt.binname));
    if opt.notbin0
        binInds(1)=[];
    end
    fieldNAll=fieldNAll(binInds);
    
    %find cells to remove for each bin
    inds2Compare=cell(1,length(fieldNAll));
    for bInd=1:length(fieldNAll)
        inds2Compare{bInd}=find(isnan(ratioDat.(fieldNAll{bInd})(:,2))); % checking only the second column is reasonable bc all col in ratioDat will be set to nan
    end
    % restrict the number of bins to consider
    inds2Compare=inds2Compare(opt.binrange);
    
    %find unique cells across all bins in inds2Compare
    catInds2Comp=[];
    for cInd=1:length(inds2Compare)
        catInds2Comp=[catInds2Comp; inds2Compare{cInd}];
    end
    catInds2Comp=unique(catInds2Comp);
    
    % now convert cell ind in all bins to nan
    for c=catInds2Comp
        for fN=1:length(fieldNAll)
            ratioDat(i).(fieldNAll{fN})(c,:)=nan;
        end %fN
    end %c
    %record number of cells removed and cells kept
    numRm=zeros(1,length(ratioDat));
    numCellTot=zeros(1,length(ratioDat));
    numRm(i)=length(catInds2Comp);
    sz=size(ratioDat(i).(fieldNAll{1}),1);
    numCellTot(i)=sz-numRm;
    
end % for i=1:lenght(ratioDat)
end% end Function