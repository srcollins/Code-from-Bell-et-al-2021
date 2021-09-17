function writeFRETMoviesTrackedStimFromImageSequences(fretImages, grayImages, filePath, pos1, varargin)
%2021-03-08 Update: This function generates a split panel movie where the left panel
%is a gray scale image and the right is a fret image. In addition, the
%function takes the pos1 shift matrix. During the experiment, cells are tracked
%and the stage coordinates for the cell centroid are saved in a variable called pos1.
%To display the cell in the lab frame of reference (where the frame is
%fixed and the cell is crawling) the image matrix needs to be shifted by the pos1
%matrix to fix the video frame. The function can also scatter the frap spot
%on the images if an optogenetic stimulation was applied. The frap
%coordinates can be added in the varargin inputs. 

%2021-03-29 Update: adding the option to flip the images left right and
%updown. The previous version could only flipup. Here are the new varargin
%calls:opt.flipud and opt.fliplr. Also added rotate90 (opt.rot90). The way that
%I configured rot90 will only let the user rotate the image once counter clockwise. 

%I also had to break the image processing and figure building into two
%loops as the counting for the frap can be i+1  depending on the
%experiment. This causes a frap spot placement error as the frap is iterated
%by i.

%Note: I have not tested wheter the flip and rotate works together.
%Importantly, the rotation is happening before either of the flips.

% pos1 is the tracked stage movements
fprintf('file path = %s\n',filePath);

%% set up defaults
opt.framerate=10;
opt.quality=100;
opt.xrange=[1 size(fretImages{1},2)];
opt.yrange=[1 size(fretImages{1},1)];
opt.boundsfret=[];   % bounds are upper and lower bounds for intensity (for display only), e.g., [100 3000]
opt.boundsgray=[];   % bounds are upper and lower bounds for intensity (for display only), e.g., [100 3000]
opt.imframerate=3;   % actual frame rate from image acquisition in seconds
opt.timezero=0;
opt.frapspot=[];
opt.numbefore=[];
opt.flipud=false;
opt.fliplr=false;
opt.rot90=false;
opt.objective=60;
opt.scalebarlength=25;
opt.scaletexttf=true;
opt.fcounter=0; % For some early driving experiments the xy coors for stage movement were recorded before the frap was applied. This causes the pos1 adj frap coors to be off by 1.


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
frapcoors=opt.frapspot+pos1 + [1-opt.xrange(1) 1-opt.yrange(1)];



vw=VideoWriter(filePath); %'\processed\'
vw.set('FrameRate',opt.framerate);
vw.set('Quality',opt.quality);
open(vw);


figure; axis image; set(gca,'XTick',[],'YTick',[]);
for i=1:length(fretImages) %process and flip/rotate images
    
    imC1{i}=sum(double(grayImages{i}),3);
    imC1{i}=[zeros(yPad1,xPad1) zeros(yPad1,imsize0(2)) zeros(yPad1,xPad2);
             zeros(imsize0(1),xPad1) imC1{i} zeros(imsize0(1),xPad2);
             zeros(yPad2,xPad1) zeros(yPad2,imsize0(2)) zeros(yPad2,xPad2)];
    
    imC1{i}=shiftMatrix(imC1{i},pos1(i,1),pos1(i,2));
    imR1{i}=fretImages{i};
    imR1{i}=[zeros(yPad1,xPad1) zeros(yPad1,imsize0(2)) zeros(yPad1,xPad2);
             zeros(imsize0(1),xPad1) imR1{i} zeros(imsize0(1),xPad2);
             zeros(yPad2,xPad1) zeros(yPad2,imsize0(2)) zeros(yPad2,xPad2)];
    imR1{i}=shiftMatrix(imR1{i},pos1(i,1),pos1(i,2));
    imR1{i}(imR1{i}==0)=nan;

    % builds the fluorescence grayscale image
    imGray=double(imC1{i}(yRangeFinal,xRangeFinal)-opt.boundsgray(1))/(opt.boundsgray(2)-opt.boundsgray(1)); %image brightness and contrast ie changing max and min on intensity histogram
    imGray=max(0,min(1,imGray));
    imGray=repmat(imGray,[1 1 3]);
    imGrayFinal{i}=imGray;
    
    imRatio=imR1{i}(yRangeFinal,xRangeFinal);
    imRatio(imRatio<opt.boundsfret(1)+.01)=opt.boundsfret(1)+.01;
    imRatio(isnan(imRatio))=opt.boundsfret(1);
    imRatio1=grs2rgb(mat2gray(imRatio,opt.boundsfret),parula_black(256));
    imRatioFinal{i}=imRatio1;
    
if opt.rot90
    imGrayFinal{i}=rot90(imGrayFinal{i});
    imRatioFinal{i}=rot90(imRatioFinal{i});
    if ~isempty(frapcoors)
        frapcoors(i,:)=fliplr(frapcoors(i,:));
        frapcoors(i,2)= (size(imRatio1,2)+1)-frapcoors(i,2);
    end
end
    
    
 

    if opt.flipud
        imDisp=flipud(imDisp);
        if ~isempty(frapcoors)
            frapcoors(i,2)=(size(imRatio,1)+1)-frapcoors(i,2);
        end
    end
    
    if opt.fliplr
        imDisp=fliplr(imDisp);
        if ~isempty(frapcoors)
            frapcoors(i,1)=(size(imRatio,2)+1)-frapcoors(i,1);
        end
    end
    
end %for i=1:length(fretImages)

for i=1:length(fretImages) %build the images
       imDisp=  [imGrayFinal{i} ones(size(imGrayFinal{i},1),2,3) imRatioFinal{i}]; %add oommented region to make grayscale and fret movies
    
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
        desiredLengthOfScale=opt.scalebarlength; % in um
        numPixInScaleBar=desiredLengthOfScale/distanceScale(opt.objective); % convert scale length in microns to pixels
        scaleHeight=opt.scalebarlength/5;
        scalePos=[0.98*figX-numPixInScaleBar, 0.85*figY];
        %place scale bar on image
        rectangle('Position',[scalePos(1) scalePos(2) numPixInScaleBar scaleHeight],'FaceColor','w','EdgeColor','w',...
            'LineWidth',3);
        %place scale bar text on image
        if opt.scaletexttf
        text(scalePos(1)+((numPixInScaleBar+5)/4), scalePos(2)+0.08*figY,sprintf('25 \x3bcm'),'Color',[1 1 1],'FontSize',14);
        end
        hold off
    end
    
    %add color bar
     
    %frap spot scatter
    if i>=opt.numbefore && ~(isempty(opt.frapspot))
        hold on;
        scatter(frapcoors(i+opt.fcounter,1)+size(imRatioFinal{i+opt.fcounter},2)+2, frapcoors(i+opt.fcounter,2),75,[1 0 1],'LineWidth',2);
        scatter(frapcoors(i+opt.fcounter,1), frapcoors(i+opt.fcounter,2),75,[1 0 1],'LineWidth',2);
    end
    %time stamp
    text(10,15,['t = ' num2str((i-opt.numbefore)*opt.imframerate) ' sec'],'Color',[1 1 1],'FontSize',12); %GB modified
    
    writeVideo(vw,getframe(gcf));
end
close(vw);
close all;