function showImagesMergeChannels(r,g,b,adjustFlag)

if nargin<4
    adjustFlag=1;
end
if nargin<3
    b=zeros(size(r));
end
if nargin<2
    g=zeros(size(r));
end

im=0*repmat(r,[1 1 3]);

if adjustFlag>0
    im=double(im);
    im(:,:,1)=mat2gray(double(r));
    im(:,:,2)=mat2gray(double(g));
    im(:,:,3)=mat2gray(double(b));
else
    im(:,:,1)=r;
    im(:,:,2)=g;
    im(:,:,3)=b;
end

imshow(im);
