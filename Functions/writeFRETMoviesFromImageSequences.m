function writeFRETMoviesFromImageSequences(fretImages, grayImages, filePath, varargin)
%2021-03-18 GB update. I added the varargin opt.tstktime to take in the metadata
%time for writing movie timestamps. I also added opt.tsstimtime to get the
%timestamps for the stimulations. if opt.tsstimtime is >1 the function will
% plot a square in the top right corner to indicate the stim duration.

% due to time constraints, I did not finish the stim indicator. 
%2021-03-30: Finishing the stim Indicator.

fprintf('file path = %s\n',filePath);

%% set up defaults
opt.framerate=10;
opt.quality=100;
opt.xrange=[1 size(fretImages{1},2)];
opt.yrange=[1 size(fretImages{1},1)];
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


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
if length(opt.xrange)==2
    opt.xrange=opt.xrange(1):opt.xrange(2);
end
if length(opt.yrange)==2
    opt.yrange=opt.yrange(1):opt.yrange(2);
end
%%
if isempty(opt.boundsgray)
    maxVals=cellfun(@(x) prctile(vect(x),99.5),grayImages);
    minVals=cellfun(@(x) prctile(vect(x),5),grayImages);
    opt.boundsgray=[mean(minVals,'omitnan') mean(maxVals,'omitnan')];
end


vw=VideoWriter(filePath); %'\processed\'
vw.set('FrameRate',opt.framerate);
vw.set('Quality',opt.quality);
open(vw);

figure; %colormap(gray)
for i=1:length(fretImages)
    
    % builds the fluorescence grayscale image    
    imGray=double(grayImages{i}(opt.yrange,opt.xrange)-opt.boundsgray(1))/(opt.boundsgray(2)-opt.boundsgray(1)); %image brightness and contrast ie changing max and min on intensity histogram
    imGray=max(0,min(1,imGray));
    imGray=repmat(imGray,[1 1 3]);
    
    imRatio=fretImages{i}(opt.yrange,opt.xrange); 
    imRatio(imRatio<opt.boundsfret(1)+.01)=opt.boundsfret(1)+.01;
    imRatio(isnan(imRatio))=opt.boundsfret(1);
    imRatio1=grs2rgb(mat2gray(imRatio,opt.boundsfret),parula_black(256));
    
    imDisp=  [imGray ones(size(imGray,1),2,3) imRatio1]; %add oommented region to make grayscale and fret movies
    
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
        scalePos=[0.98*figX-numPixInScaleBar, 0.85*figY];
        %place scale bar on image
        rectangle('Position',[scalePos(1) scalePos(2) numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
            'LineWidth',3);
        %place scale bar text on image
        if opt.scalebartexttf
        text(scalePos(1)+((numPixInScaleBar+5)/4), scalePos(2)+0.03*figY,sprintf('%s \x3bcm',num2str(opt.scalebarlength)),'Color',[1 1 1],'FontSize',14);
        end
       
    end
    
   %time stamp
   if ~isempty(opt.tstktime)
       if opt.min
       text(10,15,['t = ' num2str(opt.tstktime(i)) ' min'],'Color',[1 1 1],'FontSize',12);
       else
           text(10,25,['t = ' num2str(opt.tstktime(i)) ' s'],'Color','white',...
               'FontSize',15,'BackgroundColor','black');
       end
   else
       if opt.min
           text(10,15,['t = ' num2str((i-opt.numbefore)*opt.imframerate) ' min'],'Color',[1 1 1],'FontSize',12);
       else
           text(10,15,['t = ' num2str((i-opt.numbefore)*opt.imframerate) ' s'],'Color',[1 1 1],'FontSize',12);
       end
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