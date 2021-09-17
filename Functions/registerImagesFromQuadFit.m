function im2=registerImagesFromQuadFit(im1,p)

% Scaled parameters to work with image coordinates on a zero to one scale
% (For ease of understanding, we have the parameter values scaled such that
% p.pX(1)=1 means a 1 pixel translation for a 1024 x 1024 image
if isstruct(p)
    pX=p.pX/1024;
    pY=p.pY/1024;
else  % p is a matrix
    pX=p(:,1)/1024;
    pY=p(:,2)/1024;
end

% Convert to double data format for computations and interpolation
im1=double(im1);

%Pixel coordinates that will be used to make the coordinate grid for the image being
%aligned
[s1, s2, s3]=size(im1);
scaleFactor=max(s1,s2);
xValsOriginal=1:s2;
yValsOriginal=1:s1;
% xValsScaled=(0.5:s2)/s2 - 0.5;
% yValsScaled=(0.5:s1)/s1 - 0.5;
xValsScaled=((0.5:s2)-0.5*s2) / scaleFactor;
yValsScaled=((0.5:s1)-0.5*s1) / scaleFactor;

xMat=getXMatGridForImageRegistrationQuadFit(xValsScaled,yValsScaled);

fitdX=reshape(xMat*pX,[s1 s2]) *scaleFactor; %* s2;
fitdY=reshape(xMat*pY,[s1 s2]) *scaleFactor; %* s1;

[xCoor1,yCoor1]=meshgrid(xValsOriginal,yValsOriginal);
xCoor2=xCoor1-fitdX;   %This should give coordinates for image2 in the frame of image1
yCoor2=yCoor1-fitdY;
xCoor1rev=xCoor1+fitdX;  %This should give coordinates for image1 in the frame of image2
yCoor1rev=yCoor1+fitdY;

%Find the range of pixels in image1 that will be in the output image
xRange=[max(1,max(min(xCoor2,[],2))) min(size(im1,2),min(max(xCoor2,[],2)))];
xRange(1)=ceil(xRange(1));
xRange(2)=floor(xRange(2));
yRange=[max(1,max(min(yCoor2))) min(size(im1,1),min(max(yCoor2)))];
yRange(1)=ceil(yRange(1));
yRange(2)=floor(yRange(2));

%Find the output coordinates for the mapped image2, but in the frame of image2
xOut=xCoor1rev(yRange(1):yRange(2),xRange(1):xRange(2));
yOut=yCoor1rev(yRange(1):yRange(2),xRange(1):xRange(2));

im2(:,:,1)=im1(yRange(1):yRange(2),xRange(1):xRange(2),1);
im2(:,:,2)=interp2(xCoor1,yCoor1,im1(:,:,2),xOut,yOut);
