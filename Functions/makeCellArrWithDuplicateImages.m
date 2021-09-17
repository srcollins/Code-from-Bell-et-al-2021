function frapIms=makeCellArrWithDuplicateImages(frapIm,arrayLength)
%2020-07-18 GB: since we are now using shift matrix to adjust all of our
%images based on the stage movements between frames, I need to make an
%array of frap images that also gets the shifts applied. This function will
%take the frapImage and duplicate it in a cell array that has the length
%specified in arrayLenght. This new image array will then be processed by shiftFrapCoorsByPos1

frapIms=cell(1,arrayLength);
for i=1:arrayLength
    frapIms{i}=frapIm;
end


end