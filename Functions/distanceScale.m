function umPerPix = distanceScale(objMag, varargin)
%this function reports the microns per pixel for each objective
%magnification. The function also takes the optovar state as a second
%input. The default optovar state is "no optovar"

opt.optovar=0;
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end


if ~opt.optovar
    
    umPerPix=0.4389 * 2 * (10/objMag) *1.5;
else
    umPerPix=0.4389 * 2 * (10/objMag);
end 


end %end func