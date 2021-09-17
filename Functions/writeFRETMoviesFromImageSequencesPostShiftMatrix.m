function writeFRETMoviesFromImageSequencesPostShiftMatrix(fretImages, grayImages, filePath, frapcoors, varargin)
%% 2020-07-25 GB: This movie writing function takes images that are already
%shifted to account for stage movement between frames. Additionally, the
%function will accept a list of frap coordinates. Note, If you do not use 
%shiftMatrix data, this function will make the movie in the cell frame of
%reference where the movie pans with the cell as it crawls. This is the primary movie
%writing function for the Center Stim exp.




fprintf('file path = %s\n',filePath);

%% set up defaults
opt.framerate=10;
opt.quality=100;
opt.xrange=[1, size(fretImages{1},2)];
opt.yrange=[1, size(fretImages{1},1)];
opt.boundsfret=[];   % bounds are upper and lower bounds for intensity (for display only), e.g., [100 3000]
opt.boundsgray=[];   % bounds are upper and lower bounds for intensity (for display only), e.g., [100 3000]
opt.imframerate=3;   % actual frame rate from image acquisition in seconds
opt.timezero=0;
opt.numbefore=[];
opt.flip=false;
opt.objective=[];
opt.scalebarlength=25;


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
%%
if isempty(opt.boundsgray)
    maxVals=cellfun(@(x) prctile(vect(x),99.9),grayImages);
    minVals=cellfun(@(x) prctile(vect(x),5),grayImages);
    opt.boundsgray=[nanmean(minVals) nanmean(maxVals)];
end
if isempty(opt.boundsfret)
    maxVals=cellfun(@(x) prctile(vect(x),99),fretImages);
    minVals=cellfun(@(x) prctile(vect(x),5),fretImages);
    opt.boundsfret=[nanmean(minVals) nanmean(maxVals)];
end
% compute mapped frap coordinates

vw=VideoWriter(filePath); %'\processed\'
vw.set('FrameRate',opt.framerate);
vw.set('Quality',opt.quality);
open(vw);

figure; axis image; set(gca,'XTick',[],'YTick',[]);
for i=1:length(fretImages)
    
    imC1{i}=grayImages{i};
    
    imR1{i}=fretImages{i};
    imR1{i}(imR1{i}==0)=nan;

    % builds the fluorescence grayscale image    
    imGray=double(imC1{i}(yRangeFinal,xRangeFinal)-opt.boundsgray(1))/(opt.boundsgray(2)-opt.boundsgray(1)); %image brightness and contrast ie changing max and min on intensity histogram
    %imGray=max(0,min(1,imGray)); % this is a seanism. Not really sure how
    %it works, but it was setting the bg to 1 instead of 0, making the
    %background white
    imGray(isnan(imGray))=0;
    imGray=repmat(imGray,[1 1 3]);
    
    imRatio=imR1{i}(yRangeFinal,xRangeFinal); 
    imRatio(imRatio<opt.boundsfret(1)+.01)=opt.boundsfret(1)+.01;
    imRatio(isnan(imRatio))=opt.boundsfret(1);
    imRatio1=grs2rgb(mat2gray(imRatio,opt.boundsfret),parula_black(256));
         
    imDisp=  [imGray ones(size(imGray,1),2,3) imRatio1]; %add oommented region to make grayscale and fret movies
    if opt.flip
         imDisp=flipud(imDisp);
         frapcoors(i,2)=(size(imRatio1,1)+1)-frapcoors(i,2);
    end
    
    tickNum=linspace(0,1,4);
    tickLabels=arrayfun(@(x) {sprintf('%.2f',x)},linspace(opt.boundsfret(1),opt.boundsfret(2),4));
    h=subplot(1,1,1);
    imshow(imDisp);
    ogSize=get(gca, 'Position');
    colormap('parula'); cBar=colorbar;
    cBar.Ticks=tickNum; cBar.TickLabels=tickLabels;
    set(h, 'Position', ogSize);
    
    
    
     hold on;
    if ~isempty(opt.objective)
        [figY, figX, ~]=size(imDisp);
        distanceScale = 0.4389 * 2 * (10/opt.objective) *1.5;  %=0.21 um/pixel - converting between pixels and microns
        desiredLengthOfScale=25; % in um
        numPixInScaleBar=desiredLengthOfScale/distanceScale; % convert scale length in microns to pixels
        scaleHeight=opt.scalebarlength/5;
        scalePos=[0.98*figX-numPixInScaleBar, 0.85*figY];
        %place scale bar on image
        rectangle('Position',[scalePos(1) scalePos(2) numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
            'LineWidth',3);
        %place scale bar text on image
        text(scalePos(1)+((numPixInScaleBar+5)/4), scalePos(2)+0.03*figY,sprintf('25 \x3bcm'),'Color',[1 1 1],'FontSize',14);
        hold off
    end
    
    %add color bar
     
    %frap spot scatter
    if i>opt.numbefore 
        hold on;
        scatter(frapcoors(i,1)+size(imRatio1,2)+2, frapcoors(i,2),75,[1 0 1],'LineWidth',2);
        scatter(frapcoors(i,1), frapcoors(i,2),75,[1 0 1],'LineWidth',2);
    end
    %time stamp
    text(10,15,['t = ' num2str((i-opt.numbefore)*opt.imframerate) ' sec'],'Color',[1 1 1],'FontSize',12); %GB modified
    
    writeVideo(vw,getframe(gcf));
end
close(vw);
close all;