function writeFRETMoviesFromImageSequencesMultiPanel(mDat, filePath, varargin)
%2021-03-18 GB update. I added the varargin opt.tstktime to take in the metadata
%time for writing movie timestamps. I also added opt.tsstimtime to get the
%timestamps for the stimulations. if opt.tsstimtime is >1 the function will
% plot a square in the top right corner to indicate the stim duration.

% due to time constraints, I did not finish the stim indicator. 
%2021-03-30: Finishing the stim Indicator.

%2021-04-08: Converted to writeFRETMoviesFromImageSequencesMultiPanel. This
%function will be used to generate multipanel movies for my paper. It will
%build a movie where the gray images are postitioned above the Fret images.
%images from each folder will be positioned on the x-axis

fprintf('file path = %s\n',filePath);

%% set up defaults
opt.framerate=10;
opt.quality=100;
opt.xrange=[1 size(mDat(1).imFRET{1},2)];
opt.yrange=[1 size(mDat(1).imFRET{1},1)];
opt.boundsfret=[];   % bounds are upper and lower bounds for intensity (for display only), e.g., [100 3000]
opt.boundsgray=[];   % bounds are upper and lower bounds for intensity (for display only), e.g., [100 3000]
opt.imframerate=3;
opt.numbefore=0;
opt.objective=[];
opt.scalebarlength=25;
opt.min=false;
opt.tstktime=[];
opt.tsstimtime=[];
opt.scalebartexttf=false; % set to true if you want the text for the scale bar.
opt.labels={};

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
if length(opt.xrange)==2
    opt.xrange=opt.xrange(1):opt.xrange(2);
end
if length(opt.yrange)==2
    opt.yrange=opt.yrange(1):opt.yrange(2);
end
%% process images for movies


% check that movies are the same length - add this later. under a time
% crunch
% frameL=arrayfun(@(x) numel(x.imFRET),mDat);
% if 

%preallocate fret and grayscale cell array for final movie processing
imDim=[size(mDat(1).imFRET,1), size(mDat(1).imFRET,2), numel(mDat)];

imGrayAll=cell(imDim);
imRatioAll=cell(imDim);
tsstimtime=cell(1,imDim(3));
tstktime=cell(1,imDim(3));


for d=1:length(mDat)
    fretImages={};
    fretImages=mDat(d).imFRET;
    grayImages={};
    grayImages=mDat(d).imGray;
    
    if isempty(opt.boundsgray)
    maxVals=cellfun(@(x) prctile(vect(x),99.5),grayImages);
    minVals=cellfun(@(x) prctile(vect(x),5),grayImages);
    opt.boundsgray=[mean(minVals,'omitnan') mean(maxVals,'omitnan')];
    end
    
    tsstimtime{d}=mDat(d).time.frapTimeAdj;
    tstktime{d}=mDat(d).time.tkTimeAdj;

    
    for i=1:length(fretImages)
        
        % builds the fluorescence grayscale image
        imGray=double(grayImages{i}(opt.yrange,opt.xrange)-opt.boundsgray(1))/(opt.boundsgray(2)-opt.boundsgray(1)); %image brightness and contrast ie changing max and min on intensity histogram
        imGray=max(0,min(1,imGray));
        imGray=repmat(imGray,[1 1 3]);
        imGrayAll{1,i,d}=imGray;
        
        imRatio=fretImages{i}(opt.yrange,opt.xrange);
        imRatio(imRatio<opt.boundsfret(1)+.01)=opt.boundsfret(1)+.01;
        imRatio(isnan(imRatio))=opt.boundsfret(1);
        imRatio1=grs2rgb(mat2gray(imRatio,opt.boundsfret),parula_black(256));
        imRatioAll{1,i,d}=imRatio1;
    end
end
%% make movies

vw=VideoWriter(filePath); %'\processed\'
vw.set('FrameRate',opt.framerate);
vw.set('Quality',opt.quality);
open(vw);


vertBuffer=ones(size(imGrayAll{1},1),10,3);
horzDim=size(imGrayAll{1},2)*3+20;
horzBuffer=ones(5,horzDim,3);

figure; %colormap(gray)
  for i=1:length(fretImages)  
      
     %Build composite image 
    imDisp=  [imGrayAll{1,i,1}, vertBuffer imGrayAll{1,i,2}, vertBuffer imGrayAll{1,i,3};...
        horzBuffer;...
        imRatioAll{1,i,1} vertBuffer imRatioAll{1,i,2} vertBuffer imRatioAll{1,i,3}]; %add oommented region to make grayscale and fret movies
    
    % set colorbar
    tickNum=linspace(0,1,4);
    tickLabels=arrayfun(@(x) {sprintf('%.2f',x)},linspace(opt.boundsfret(1),opt.boundsfret(2),4));
    h=subplot(1,1,1);
    imshow(imDisp);
    ogSize=get(gca, 'Position');
    colormap('parula'); cBar=colorbar;
    cBar.Ticks=tickNum; cBar.TickLabels=tickLabels;
    set(h, 'Position', ogSize);
    
    %set scalebar
     hold on;
    if ~isempty(opt.objective)
        [figY, figX, ~]=size(imDisp);
        umPerPix=distanceScale(opt.objective);  %=0.21 um/pixel - converting between pixels and microns
        desiredLengthOfScale=opt.scalebarlength; % in um
        numPixInScaleBar=desiredLengthOfScale/umPerPix; % convert scale length in microns to pixels
        scaleHeight=opt.scalebarlength/5;
        scalePos=[0.98*figX-numPixInScaleBar, 0.95*figY];
        %place scale bar on image
        rectangle('Position',[scalePos(1) scalePos(2) numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
            'LineWidth',3);
        %place scale bar text on image
        if opt.scalebartexttf
        text(scalePos(1)+((numPixInScaleBar+5)/4), scalePos(2)+0.03*figY,sprintf('%s \x3bcm',num2str(opt.scalebarlength)),'Color',[1 1 1],'FontSize',14);
        end
       
    end
    
   %time stamp 
   spaceDist=size(imGrayAll{1,i,1},2) + size(vertBuffer,2);
   
   for m=1:length(tstktime)
       xPos=10+(m-1)*spaceDist;
       textIn=sprintf('t= %.1f s',tstktime{m}(i));
       text(xPos,25,textIn,'Color','white',...
           'FontSize',10,'BackgroundColor','black');
       
       if ~isempty(tsstimtime{m})
           if tstktime{m}(i)>tsstimtime{m}(1) && tstktime{m}(i)<tsstimtime{m}(end)
               stimPostition=[xPos,45,25,25];
               rectangle('Position',stimPostition,'FaceColor','w','EdgeColor','w');
               
           end
       end
       
       text(xPos+250,25,opt.labels{m},'Color','white',...
           'FontSize',10,'BackgroundColor','black');
   end
   
   
   
%    % stim indicator for multistim pdl exp
   if ~isempty(opt.tsstimtime)
       if opt.tstktime(i)>opt.tsstimtime(1) && opt.tstktime(i)<opt.tsstimtime(end)
           stimPostition=[10,45,25,25];
           rectangle('Position',stimPostition,'FaceColor','w','EdgeColor','w');
           
       end
   end
   hold off;
   writeVideo(vw,getframe(gcf));
end
close(vw);
close all;