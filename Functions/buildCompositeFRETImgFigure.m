function buildCompositeFRETImgFigure(fretImages, filePath,frameSelect, varargin)
%% This function will generate a composite fret image figure. And optionally, a gray image figure
% The function will take a series of fret images and build one composite
% image that is a montage of the input images. The function will also 
% scatter the frap spot and display the color bar.
% 
% Inputs: The function requires a cell array of fret images, the filepath,
% and the vector containing the two image frames.  
%The function will also give the user a number of name/value pair calls,
% including the option to include pos1 for shift matrix and the ability to
%add the gray scale images to the final figure. See the additional
% options below. 
%updated 2020-04-07 GB

%updated 2021-01-30: fixed scale bar bug, also added an option to create a
%ratio diff figure. This will be saved as a separate pdf

% 2021-02-03: added the option to create a ratio diff image that is stored
% in a separate pdf.

% 2021-02-04: add optional position and units to imshow so that the user
% can specify how large the figure is. 

%2021-03-06: added flipud and fliplr options.

%2021-03-18: added the option to not show the colorbars
%2021-04-06: allow differnce image output to be relative fret ratio
%FrameAfter./FrameBefore; use opt.relchangetf
fprintf('file path = %s\n',filePath);

%% set up defaults
opt.xrange=[1 size(fretImages{1},2)];
opt.yrange=[1 size(fretImages{1},1)];
opt.boundsfret=[];   % bounds are upper and lower bounds for intensity (for display only), e.g., [100 3000]
opt.imframerate=3;   % actual frame rate from image acquisition in seconds
opt.ts=[];
opt.timezero=0;
opt.frapspot=[];
opt.numbefore=[];
opt.flipud=false;
opt.fliplr=false;
opt.objective=[];
opt.scalebarlength=25;
opt.savefiletype='tif';
opt.diffrange=[-0.2 0.2]; % set the bounds for the diff image
opt.pos1=[];
opt.imsgray={};
opt.boundsgray=[];
opt.diffimgtf=false;
opt.difframes=[1,4];
opt.diffcrange=[-.1 .1];
opt.frapscattersize=75;
opt.scalebartexttf=true;
opt.colorbartf=true;
opt.relchangetf=false; % use relative chagne in fret (FRET-imAfter./FretImgBefore)
opt.rcrange=[0.99 1.09]; % relative change in fret signaling range.

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
if length(opt.xrange)==2
    opt.xrange=opt.xrange(1):opt.xrange(2);
end
if length(opt.yrange)==2
    opt.yrange=opt.yrange(1):opt.yrange(2);
end

imsize0=size(fretImages{1});
xPad1=max(0,1-opt.xrange(1));
xPad2=max(0,opt.xrange(end)-imsize0(2));
xRangeFinal=(opt.xrange)+xPad1;
yPad1=max(0,1-opt.yrange(1));
yPad2=max(0,opt.yrange(end)-imsize0(1));
yRangeFinal=(opt.yrange)+yPad1;

imRatios=cell(size(frameSelect));
imRatios1=cell(size(frameSelect));
%%
if isempty(opt.boundsgray)
    maxVals=cellfun(@(x) prctile(vect(x),99.9),opt.imsgray);
    minVals=cellfun(@(x) prctile(vect(x),5),opt.imsgray);
    opt.boundsgray=[nanmean(minVals) nanmean(maxVals)];
end

if isempty(opt.boundsfret)
    maxVals=cellfun(@(x) prctile(vect(x),99),fretImages);
    minVals=cellfun(@(x) prctile(vect(x),5),fretImages);
    opt.boundsfret=[nanmean(minVals) nanmean(maxVals)];
end
% compute mapped frap coordinates
if ~isempty(opt.pos1)
frapcoors=opt.frapspot+opt.pos1 + [1-opt.xrange(1) 1-opt.yrange(1)];
elseif isempty(opt.frapspot)
    frapcoors=[];
else
frapcoors=opt.frapspot + [1-opt.xrange(1) 1-opt.yrange(1)];
end

%build fret ratio images
count=0;
for i=frameSelect
    count=count+1;
    
    imR1=fretImages{i};
    imR1=[zeros(yPad1,xPad1) zeros(yPad1,imsize0(2)) zeros(yPad1,xPad2);
        zeros(imsize0(1),xPad1) imR1 zeros(imsize0(1),xPad2);
        zeros(yPad2,xPad1) zeros(yPad2,imsize0(2)) zeros(yPad2,xPad2)];
    if ~isempty(opt.pos1)
    imR1=shiftMatrix(imR1,opt.pos1(i,1),opt.pos1(i,2));
    end
    imR1(imR1==0)=nan;

    imRatio=imR1(yRangeFinal,xRangeFinal); 
    if opt.flipud
        imRatio=flipud(imRatio);
        if ~isempty(frapcoors)
            frapcoors(i,2)=(size(imRatio,1)+1)-frapcoors(i,2);
        end
    end
    
    if opt.fliplr
        imRatio=fliplr(imRatio);
        if ~isempty(frapcoors)
            frapcoors(i,1)=(size(imRatio,2)+1)-frapcoors(i,1);
        end
    end
    
    imRatio4Diff{count}=imRatio;
    imRatio(imRatio<opt.boundsfret(1)+.01)=opt.boundsfret(1)+.01;
    imRatio1=imRatio;
    imRatio1(isnan(imRatio1))=opt.boundsfret(1);
    imRatios{count}=imRatio;
    imRatios1{count}=imRatio1;
end % i=frameSelect

% build Ratio Diff image
if opt.diffimgtf || opt.relchangetf
    if opt.diffimgtf
    diffImg=imRatio4Diff{opt.difframes(2)}-imRatio4Diff{opt.difframes(1)};
    diffImg(diffImg<opt.diffcrange(1)+.001)=opt.diffcrange(1)+.001;
    diffImg1=diffImg;
    diffImg1(isnan(diffImg1))=opt.diffcrange(1);
    diffImg2=grs2rgb(mat2gray(diffImg1,opt.diffcrange),blue_yellow_gray(256));
    imDispDiff=diffImg2;
    elseif opt.relchangetf % relative change instead of difference
        diffImg=imRatio4Diff{opt.difframes(2)}./imRatio4Diff{opt.difframes(1)};
    diffImg(diffImg<opt.rcrange(1)+.01)=opt.rcrange(1)+.01;
    diffImg1=diffImg;
    diffImg1(isnan(diffImg1))=opt.rcrange(1);
    diffImg2=grs2rgb(mat2gray(diffImg1,opt.rcrange),blue_yellow_gray(256));
    imDispDiff=diffImg2;
    end
    
    if opt.flipud
         imDisp=flipud(imDispDiff);
         if ~isempty(frapcoors)
         frapcoors(i,2)=(size(imRatio1,1)+1)-frapcoors(i,2);
         end
    end
    
    if opt.fliplr 
         imDisp=fliplr(imDispDiff);
         if ~isempty(frapcoors)
         frapcoors(i,1)=(size(imRatio1,2)+1)-frapcoors(i,1);
         end
    end
    
    tickNum=linspace(0,1,4);
    if opt.diffimgtf
        tickLabels=arrayfun(@(x) {sprintf('%.2f',x)},linspace(opt.diffcrange(1),opt.diffcrange(2),4));
    elseif opt.relchangetf
        tickLabels=arrayfun(@(x) {sprintf('%.2f',x)},linspace(opt.rcrange(1),opt.rcrange(2),4));
    end
        
    h=subplot(1,1,1);
    imshow(imDispDiff);
    if opt.colorbartf
    ogSize=get(gca, 'Position');
    colormap('blue_yellow_gray'); cBar=colorbar;
    cBar.Ticks=tickNum; cBar.TickLabels=tickLabels;
    set(h, 'Position', ogSize);
    end
    hold on;
    %scale bar
    [figY, figX, ~]=size(imDispDiff);
    umPerPix=distanceScale(opt.objective);  %=0.21 um/pixel - converting between pixels and microns
    numPixInScaleBar=opt.scalebarlength/umPerPix; % convert scale length in microns to pixels
    scaleHeight=opt.scalebarlength/5;
    scalePos=[0.98*figX-numPixInScaleBar, 0.85*figY];
    %place scale bar on image
    rectangle('Position',[scalePos(1) scalePos(2) numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
        'LineWidth',3);
    %place scale bar text on image
    if opt.scalebartexttf
    text(scalePos(1)+((numPixInScaleBar+5)/4), scalePos(2)+0.03*figY,sprintf('%s \x3bcm', num2str(opt.scalebarlength)),'Color',[1 1 1],'FontSize',14);
    end
    %Frap
    scatter(frapcoors(opt.difframes(1),1), frapcoors(opt.difframes(1),2),opt.frapscattersize,[1 0 1],'LineWidth',2);
    hold off;
    if opt.diffimgtf
    ratioDiffName=sprintf('ratioDiff%s-%s.pdf',num2str(frameSelect(opt.difframes(2))),num2str(frameSelect(opt.difframes(1))));
    elseif opt.relchangetf
        ratioDiffName=sprintf('relChangeFret%s-%s.pdf',num2str(frameSelect(opt.difframes(2))),num2str(frameSelect(opt.difframes(1))));
    end
    saveRootWithFileName=fullfile(filePath,ratioDiffName);
    export_fig(saveRootWithFileName)
    
end

% build Gray images
if ~isempty(opt.imsgray)
    count=0;
    for i=frameSelect
        count=count+1;
        
        imC1=opt.imsgray{i};
        imC1=[zeros(yPad1,xPad1) zeros(yPad1,imsize0(2)) zeros(yPad1,xPad2);
            zeros(imsize0(1),xPad1) imC1 zeros(imsize0(1),xPad2);
            zeros(yPad2,xPad1) zeros(yPad2,imsize0(2)) zeros(yPad2,xPad2)];
        if ~isempty(opt.pos1)
            imC1=shiftMatrix(imC1,opt.pos1(i,1),opt.pos1(i,2));
        end
        % builds the fluorescence grayscale image    
    imGray=double(imC1(yRangeFinal,xRangeFinal)-opt.boundsgray(1))/(opt.boundsgray(2)-opt.boundsgray(1)); %image brightness and contrast ie changing max and min on intensity histogram
    imGray(isnan(imGray))=-1;
    imGray=max(0,min(1,imGray));
     if opt.flipud
        imGray=flipud(imGray);
        if ~isempty(frapcoors)
            frapcoors(i,2)=(size(imRatio,1)+1)-frapcoors(i,2);
        end
    end
    
    if opt.fliplr
        imGray=fliplr(imGray);
        if ~isempty(frapcoors)
            frapcoors(i,1)=(size(imRatio,2)+1)-frapcoors(i,1);
        end
    end
    
    imGray=repmat(imGray,[1 1 3]);
    imsGray{count}=imGray;
    end %i frameselct
end %if ~isempty(opt.imsgray)

%build ratio composite for display
imRatios1=cellfun(@(x) grs2rgb(mat2gray(x,opt.boundsfret),parula_black(256)),imRatios1,'UniformOutput',false);
imgGap=ones(size(imRatios1{1},1),2,3);
for imD=1:length(frameSelect)
    if imD==1
imDisp= [imRatios1{imD}]; %[imRatios1{imD},imgGap];
    else
        imDisp=[imDisp, imgGap imRatios1{imD}];
    end
end
   
%build grayscale composite img for display
if ~isempty(opt.imsgray)
    for imG=1:length(frameSelect)
        if imG==1
            imDispGray= [imsGray{imG}]; %[imsGray{imG},imgGap]
        else
            imDispGray=[imDispGray, imgGap imsGray{imG}];
        end
    end
end



% if opt.flipud
%     imDisp=flipud(imDisp);
%     imDispGray=flipud(imDispGray);
%     if ~isempty(frapcoors)
%         frapcoors(i,2)=(size(imRatio1,1)+1)-frapcoors(i,2);
%     end
% end
% 
% if opt.fliplr
%     imDisp=fliplr(imDisp);
%     imDispGray=fliplr(imDispGray);
%     if ~isempty(frapcoors)
%         frapcoors(i,1)=(size(imRatio1,2)+1)-frapcoors(i,1);
%     end
% end

% show ratio image composite and add timer, scalebar frapspot and colorbar
tickNum=linspace(0,1,4);
tickLabels=arrayfun(@(x) {sprintf('%.2f',x)},linspace(opt.boundsfret(1),opt.boundsfret(2),4));
h=subplot(1,1,1);
imshow(imDisp);
 if opt.colorbartf
ogSize=get(gca, 'Position');
colormap('parula_black'); cBar=colorbar;
cBar.Ticks=tickNum; cBar.TickLabels=tickLabels;
set(h, 'Position', ogSize);
 end




hold on;
if ~isempty(opt.objective)
    [figY, figX, ~]=size(imDisp);
    umPerPix=distanceScale(opt.objective);  %=0.21 um/pixel - converting between pixels and microns
    numPixInScaleBar=opt.scalebarlength/umPerPix; % convert scale length in microns to pixels
    scaleHeight=opt.scalebarlength/5;
    scalePos=[0.98*figX-numPixInScaleBar, 0.85*figY];
    %place scale bar on image
    rectangle('Position',[scalePos(1) scalePos(2) numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
        'LineWidth',3);
    %place scale bar text on image
    if opt.scalebartexttf
    text(scalePos(1)+((numPixInScaleBar+5)/4), scalePos(2)+0.03*figY,sprintf('%s \x3bcm', num2str(opt.scalebarlength)),'Color',[1 1 1],'FontSize',14);
    end
    hold off
end



%frap spot scatter
imgSpacer=size(imRatios1{1},2);

    for imtext=1:length(frameSelect)
        hold on;
        %add timer
        multiplier=[1:length(frameSelect)]-1;
        if isempty(opt.ts)
            timeStamp=num2str((frameSelect(imtext)-opt.numbefore)*opt.imframerate);
        else
            timeStamp=num2str(round(opt.ts.timesInfrapRef(frameSelect(imtext)),1));
        end
        textVals=sprintf('t = %s sec',timeStamp);
        text(10+((size(imgGap,2)+imgSpacer)*multiplier(imtext)),15,...
            textVals,'Color',[1 1 1],'FontSize',12); % gray1 timestap
        %add frapspot
        if frameSelect(imtext)>opt.numbefore && ~isempty(opt.frapspot)
            hold on
            scatter(frapcoors(frameSelect(imtext)+1,1)+((size(imgGap,2)+imgSpacer)*multiplier(imtext)), frapcoors(frameSelect(imtext)+1,2),opt.frapscattersize,[1 0 1],'LineWidth',2);
        end
    end
   
    
    if strcmpi(opt.savefiletype,'tif')
        imRatioName=sprintf('FretCompositeframe%s.tif',num2str(frameSelect));
        saveRootWithFileName=fullfile(filePath,imRatioName);
        export_fig(saveRootWithFileName)
    elseif strcmpi(opt.savefiletype, 'pdf')
        imRatioName=sprintf('FretCompositeframe%s.pdf',num2str(frameSelect));
        saveRootWithFileName=fullfile(filePath,imRatioName);
        export_fig(saveRootWithFileName)
        %printFigurePDF(saveRootWithFileName,'-bestfit');
    end



%make gray scale image ---------------------------
% show gray image composite and add timer, scalebar frapspot and colorbar
if ~isempty(opt.imsgray)
tickNum=linspace(0,1,4);
tickLabels=arrayfun(@(x) {sprintf('%.2f',x)},linspace(floor(opt.boundsgray(1)),floor(opt.boundsgray(2)),4));
h=subplot(1,1,1);
imshow(imDispGray);
 if opt.colorbartf
ogSize=get(gca, 'Position');
colormap('gray'); cBar=colorbar;
cBar.Ticks=tickNum; cBar.TickLabels=tickLabels;
set(h, 'Position', ogSize);
 end


hold on;
if ~isempty(opt.objective)
    [figY, figX, ~]=size(imDisp);
    umPerPix=distanceScale(opt.objective);  %=0.21 um/pixel - converting between pixels and microns
    numPixInScaleBar=opt.scalebarlength/umPerPix; % convert scale length in microns to pixels
    scaleHeight=opt.scalebarlength/5;
    scalePos=[0.98*figX-numPixInScaleBar, 0.85*figY];
    %place scale bar on image
    rectangle('Position',[scalePos(1) scalePos(2) numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
        'LineWidth',3);
    %place scale bar text on image
    if opt.scalebartexttf
    text(scalePos(1)+((numPixInScaleBar+5)/4), scalePos(2)+0.03*figY,sprintf('%s \x3bcm',num2str(opt.scalebarlength)),'Color',[1 1 1],'FontSize',14);
    end
    hold off
end



%frap spot scatter
imgSpacer=size(imRatios1{1},2);

   for imtext=1:length(frameSelect)
        hold on;
        %add timer
        multiplier=[1:length(frameSelect)]-1;
        if isempty(opt.ts)
        timeStamp=num2str((frameSelect(imtext)-opt.numbefore)*opt.imframerate);
        else
            timeStamp=num2str(round(opt.ts.timesInfrapRef(frameSelect(imtext)),1));
        end
        textVals=sprintf('t = %s sec',timeStamp);
        text(10+((size(imgGap,2)+imgSpacer)*multiplier(imtext)),15,...
            textVals,'Color',[1 1 1],'FontSize',12); % gray1 timestap
        %add frapspot
        if frameSelect(imtext)>opt.numbefore && ~isempty(opt.frapspot)
            hold on
            scatter(frapcoors(frameSelect(imtext)+1,1)+((size(imgGap,2)+imgSpacer)*multiplier(imtext)), frapcoors(frameSelect(imtext)+1,2),opt.frapscattersize,[1 0 1],'LineWidth',2);
        end
    end
    
    
    if strcmpi(opt.savefiletype,'tif')
        imRatioName=sprintf('grayComposite%s.tif',num2str(frameSelect));
        saveRootWithFileName=fullfile(filePath,imRatioName);
        export_fig(saveRootWithFileName)
    elseif strcmpi(opt.savefiletype, 'pdf')
        imRatioName=sprintf('grayComposite%s.pdf',num2str(frameSelect));
        saveRootWithFileName=fullfile(filePath,imRatioName);
        export_fig(saveRootWithFileName)
        %printFigurePDF(saveRootWithFileName,'-bestfit');
    end
end


%close all;