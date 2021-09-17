function errVals=getErrorValuesForRegisteredImages(p,imStack)

imStack=registerImagesFromQuadFit(imStack,p);
s=size(imStack);                    % The size of the output image is used to 
                                    % scale the error values so that there is 
                                    % no computational advantage of generating a smaller image
% factor=sqrt(s(1)*s(2));
% errVals=vect(imStack(:,:,1)-imStack(:,:,2))/factor;
% errVals=sum(vect(errVals.*errVals));
errVals=1-myNanCorrcoef(vect(imStack(:,:,1)),vect(imStack(:,:,2)));
