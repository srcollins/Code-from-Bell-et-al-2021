function imOut = imageBin(im,binsize)
%IMAGEBIN uses binning to make a smaller image. It throws away border
%pixels that don't land in a full bin.

[nr nc n3]=size(im);
nrOut=floor(nr/binsize);
ncOut=floor(nc/binsize);

im=double(im);
imOut=zeros(nrOut,ncOut,n3);
for i=1:binsize
    for j=1:binsize
        imOut=imOut+im(i:binsize:(nrOut*binsize),j:binsize:(ncOut*binsize),:);
    end
end
imOut=imOut/binsize/binsize;

end

