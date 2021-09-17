function imgsOut=convertZeroPadsToNans(imgIn,masks)
%2020-07-25 GB: This function is for images that were just shifted using shiftCellImWithPos1
%to account for stage movement. Masks is assumed to already be shifted with
%shiftCellImWithPos1. This function will use masks to id pixels in imgIn
%that are also in the mask, and will set all other pixes to nan. If imgIn
%was shifted this function will convert the zeros that were used to pad
%imgIn with nans.
imgsOut=cell(1,length(imgIn));

for i=1:length(imgIn)
    tempIm=imgIn{i};
    tempIm(~masks{i})=nan;
    imgsOut{i}=tempIm;
end