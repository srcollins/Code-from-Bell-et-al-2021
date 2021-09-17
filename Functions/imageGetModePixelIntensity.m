function mode=imageGetModePixelIntensity(im,bins,samplingRatio)
% Computes the mode pixel value for an image, using ksdensity. The second
% input "bins" can be used to specify which intensity values are examined

if nargin<2
    bins=0:5:5000;
end
if nargin<3
    samplingRatio=10;
end
im=double(im);
[f,xi]=ksdensity(im(1:samplingRatio:end),bins);
[~,ind]=max(f);
mode=xi(ind);
