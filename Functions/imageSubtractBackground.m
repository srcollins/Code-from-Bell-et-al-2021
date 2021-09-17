function im1=imageSubtractBackground(im,p_tile,blocksize)
%creates a new image background subtracted image. It uses the prctile
%percentile pixel intensity in a blocksize x blocksize neighborhood as the background
%intensity. prctile should be specified as an integer (e.g. 10 = 10th
%percentile).

if nargin<2
    p_tile=20;    %default - use the 10th percentile
end
if nargin<3
    blocksize=32;
end

if ndims(im)==3
    [r c s]=size(im);
else
    [r c]=size(im);
end
if ndims(im)==2 || s==1     %single image
    im1=processSingleImage(im,p_tile,blocksize);
end
if ndims(im)==3 && s>1      %image stack
    im1=im*0;
    for j=1:s
        im1(:,:,j)=processSingleImage(im(:,:,j),p_tile,blocksize);
    end    
end
end

%------------------------------------------------------------------------

function im1=processSingleImage(im,p_tile,blocksize)
%this subfunction does the work (one for one image)

im=squeeze(im); %remove the singleton dimension if one exists

% Remove any horizontal stripes of background -- This was an artifact I needed to handle with an older camera. It can be uncommented if needed.
% im=im-uint16(repmat(prctile(double(im),15,2),[1 size(im,2)]));

%Pad the image to make it divisible by blocksize
temp=nan(ceil(size(im,1)/blocksize)*blocksize,ceil(size(im,1)/blocksize)*blocksize);
[t1 t2]=size(temp);
[i1 i2]=size(im);
s1=ceil((t1-i1)/2);
s2=ceil((t2-i2)/2);
% temp(s1:(s1+i1-1),s2:(s2+i2-1))=im; Gives error bc first index is 0 not 1
temp((s1+1):(s1+i1),(s2+1):(s2+i2))=im;
temp(s1:-1:1,:)=temp(s1:(2*s1-1),:);
temp(end:-1:(end-s1-1),:)=temp((end-s1*2-1):(end-s1),:);
temp(:,s2:-1:1)=temp(:,s2:(2*s2-1),:);
temp(:,end:-1:(end-s2-1))=temp(:,(end-s2*2-1):(end-s2));

% Do the local computation
temp=blkproc(temp,[blocksize blocksize],@(x) prctile(x(:),p_tile));
temp=imresize(temp,[t1 t2],'bilinear');
% bg=temp(s1:(s1+i1-1),s2:(s2+i2-1));
bg=temp((s1+1):(s1+i1),(s2+1):(s2+i2));
%bg=imresize(temp,size(im),'bilinear');
im1=double(im)-bg;

end
