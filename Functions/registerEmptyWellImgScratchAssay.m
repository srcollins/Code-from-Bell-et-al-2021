function [regEmptyYFP, regEmptyCFP] = registerEmptyWellImgScratchAssay(emptyYFP,emptyCFP,p)
%% This function is desgined to receive two images, register them with the
%p alignment parameters, and them output the individual images.
allEmpty(:,:,1)=emptyYFP;
allEmpty(:,:,2)=emptyCFP;
regImgs=registerImagesFromQuadFit(allEmpty,p);
regEmptyYFP=regImgs(:,:,1);
regEmptyCFP=regImgs(:,:,2);


end