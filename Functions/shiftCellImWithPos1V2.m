function shiftedCells = shiftCellImWithPos1V2(cellImg,pos1,varargin)
%2020-07-16 GB: This function is designed to take cell arrays of cell
%Images and shift the cells to account for the stage translation stored in
%pos1. The pos1 input should be a list of x,y coordinates that are in units
%of pixels. 

% Inputs: cell array of images that are size n by m by 1

%2020-07-25 Update GB: 1) added logical as varargin input allowing the
%image to be converted from double to logical. This is used for the
%frapMask and Mask images. 

%2020-07-28 V2 update: after talking with sean, this function does not need
%the padding. Also the shifting should only occur for when the stage moved.
%Finally we will consider frame 6 the reference frame, and the shifting
%will be relative to that frame.
%% process varargin inputs
opt.numbefore=5;
opt.logical=false;

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%% Align images based on frame 6 (which should have the correct FRAP spot location)

% for i=1:length(filenames)
%     load([folder filesep filenames{i}]);
%     imDat.fretFinal=imDat.cropFiltFret;
%     imDat.centFinal=imDat.cellCentroids;
%     for j=1:5  %the other images should not need to be shifted for alignment
%         imDat.fretFinal{j}=shiftMatrix(imDat.cropFiltFret{j},imDat.pos1(j,1)-imDat.pos1(6,1),imDat.pos1(j,2)-imDat.pos1(6,2));
%         imDat.fretFinal{j}(imDat.fretFinal{j}==0)=nan;
%         imDat.centFinal(j,:)=imDat.cellCentroids(j,:)+[imDat.pos1(j,1)-imDat.pos1(6,1) imDat.pos1(j,2)-imDat.pos1(6,2)];
%     end
%     save([folder filesep filenames{i}],'imDat');    
% end

shiftedCells=cellfun(@(x) x,cellImg,'UniformOutput',false);


    % shift the images
    for i=1:opt.numbefore
        shiftedCells{i}=shiftMatrix(shiftedCells{i},pos1(i,1)-pos1(6,1),pos1(i,2)-pos1(6,2));
        
        if opt.logical
            shiftedCells{i}=logical(shiftedCells{i});
        end
    end

end
