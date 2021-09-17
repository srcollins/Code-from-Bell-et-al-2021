function p=computeRegistrationParametersQuadFit(imRef,p0)
% Using nonlinear optimization to compute registration parameters to align
% imRef(:,:,2) to imRef(:,:,1). Starting parameter guesses can optionally
% be set using p0. If p0 is not passed in to the function, it will begin
% with a translation only registration to get the translation approximately
% right with a small parameter space. If p0 is passed in, the function will
% go straight to full optimization.

% Generate Normalized Reference Image and Gradient Image
imgA=imRef;
mode1=imageGetModePixelIntensity(imRef(:,:,1),0:round(prctile(vect(imRef(:,:,1)),90)));
imgA(:,:,1)=imRef(:,:,1)-mode1;
imgA(:,:,1)=imgA(:,:,1)/mean(vect(imgA(:,:,1)));
mode2=imageGetModePixelIntensity(imRef(:,:,2),0:round(prctile(vect(imRef(:,:,2)),90)));
imgA(:,:,2)=imRef(:,:,2)-mode2;
imgA(:,:,2)=imgA(:,:,2)/mean(vect(imgA(:,:,2)));
imgG(:,:,1)=imgradient(imfilter(imgA(:,:,1),fspecial('gaussian',5,2),'symmetric'));
imgG(:,:,2)=imgradient(imfilter(imgA(:,:,2),fspecial('gaussian',5,2),'symmetric'));

% If p0 is not passed in, begin with: INITIAL TRANSLATION-ONLY REGISTRATION
if nargin<2 || isempty(p0)
    bin1=round(max(size(imgA))/256); bin1=max(bin1,1);
    imgA1=imageBin(imgA,bin1);
    f=@(x)getErrorValuesForRegisteredImages([x; zeros(5,2)],imgA1);
    p0=[0.1 0.1];
    tic;
    pFit=fminunc(f,p0); disp(pFit(end,:));
    pFit=[pFit; zeros(5,2)];
    toc;  % Display the time taken for the translation-only registration
    fprintf('Error after translation-only registration: %.4f\n', getErrorValuesForRegisteredImages(pFit,imRef)); % Display error magnitude after translation-only registration
else
    if isstruct(p0)
        p0=[p0.pX(:) p0.pY(:)];  % For computation in this function, we reorganize the parameters as a matrix. They get put back into a structure at the end of the function
    end
    fprintf('Error with initial parameter guess: %.4f\n', getErrorValuesForRegisteredImages(p0,imRef)); % Display error magnitude after translation-only registration
    pFit=p0;  % Start the next steps using the initial parameter guess
end

% INITIAL CRUDE REGISTRATION
% an initial crude registration is performed using a binned, lower
% resolution image. This should make computation faster and get us close to
% the optimal parameter values.
bin1=round(max(size(imgA))/256); bin1=max(bin1,1);
imgA1=imageBin(imgA,bin1);
f=@(x)getErrorValuesForRegisteredImages(x,imgA1);
p0=pFit;
tic;
options=optimset('MaxIter',5000,'MaxFunEvals',5000);
pFit=fminsearch(f,p0,options); disp(pFit(end,:));  % Derivative-free method
toc;
fprintf('Error after 1st crude all-parameter registration: %.4f\n', getErrorValuesForRegisteredImages(pFit,imRef)); 
tic;
options=optimset('MaxIter',5000,'MaxFunEvals',5000);
pFit=fminunc(f,pFit,options); disp(pFit(end,:));    % More standard gradient-based approach - zero in on the minimum
toc;
fprintf('Error after 2nd (gradient-based) crude all-parameter registration: %.4f\n', getErrorValuesForRegisteredImages(pFit,imRef)); % Display error magnitude after translation-only registration

% REPEAT WITH HIGHER IMAGE RESOLUTION
% Secondary more precise registration using edge detection
bin2=round(max(size(imgA))/512); bin2=max(bin2,1);
if bin2<bin1
    imgG0=imageBin(imgG,bin2);
    f=@(x)getErrorValuesForRegisteredImages(x,imgG0);
    options=optimset('MaxIter',1000,'MaxFunEvals',1000);
    tic;
    %pFit=fminunc(f,pFit,options); disp(pFit(end,:));
    pFit=fminsearch(f,pFit); disp(pFit(end,:));
    toc;
    fprintf('Error after final all-parameter registration: %.4f\n', getErrorValuesForRegisteredImages(pFit,imRef));
end
p.pX=pFit(:,1); p.pY=pFit(:,2);
