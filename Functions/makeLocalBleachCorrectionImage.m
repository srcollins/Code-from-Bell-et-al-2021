function bleachCorrectionImage=makeLocalBleachCorrectionImage(p,distMask,relativeLaserPower,t)
% This function computes a bleaching correction image that can be used to
% correct observed data.
%
% Inputs:
% p are the bleaching correction parameters
% p(1) is a magnitude parameter which should be appropriate for a reference
% laser power - initially defined as the 15 mW setting as used in July,
% 2020
% p(2) is a shape parameter which describes the initial bleaching pattern.
% distMask should be precomputed as a matrix indicating the distance in
% microns from every spot in the image to the center of the laser spot
% t is time since the laser pulse

%20201004 update GB: if the bleach correction image is all zeros resulting
%from negative time, convert the image to ones of the same size. 

p(1)=p(1)*relativeLaserPower;
bleachCorrectionImage=1-bleachedMoleculeDiffusionFunction(p,distMask,t);
if all(vect(bleachCorrectionImage)==0)
    bleachCorrectionImage=double(ones(size(bleachCorrectionImage,1),size(bleachCorrectionImage,2)));
end
