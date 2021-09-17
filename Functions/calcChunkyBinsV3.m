function cBinOut = calcChunkyBinsV3(imAcceptor,imDonor,frapMask,ratioCorrectionFinal,bleachCorrection, varargin)
% 2020-05-28 GB: This function is designed to calculate the whole cell mean
% (where the mean of each channel is calculated before the ratio is taken.
% Additionally, the function will calculate a distance mask for each image
% where the pixels are measured from the fixed frap spot location. Three
% bins will be created (0-3um,3-6um,6-9um). Pixels in each bin with be used
% to calculate the means for each channel and the mean fret ratio for each
% bin will be calculated.

%update 2020-07-23 GB: frapMask can now take a cell array of masks rather
%than a single image. This allows for pos1 adjusted frap images.

%udpdate 2020-08-05 GB: since we are using the ratioCorrectionFinal image
%after the other images have been shifted, I also need to use shifted
%ratioCorrection images. I have modified the script to handle the
%shifted ratioCorrection images as a cell array. 

%2020-09-02 V3 update: added in the bleaching correction. Also added
%stimInd. Stim Ind should be a logical vector with zeros pre stimulus and
%ones post stim.
%% process varargin inputs
opt.objmag=60;
opt.binedges=[3,6,9];
opt.numbefore=5;

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% process images for chunkyBin!
%first calc the whole cell mean fret ratios
% corrYFP=cellfun(@(x,y) x./y,imAcceptor,ratioCorrectionFinal,'UniformOutput',false); %apply the ratio correction.
% corrYFP=cellfun(@(x,y) double(x./y),corrYFP,bleachCorrection,'UniformOutput',false); %apply the bleaching correction.
% % apply bleach correction as a for loop
corrYFP=cell(1,length(imAcceptor));
for bc=1:length(bleachCorrection)
    corrYFP{bc}=imAcceptor{bc}./ratioCorrectionFinal{bc};
    corrYFP{bc}=corrYFP{bc}./bleachCorrection{bc};
end
meanYFP=cellfun(@(x) nanmean(x,'all'),corrYFP,'UniformOutput',false);
meanCFP=cellfun(@(x) nanmean(x, 'all'),imDonor,'UniformOutput',false);
wholeCellMeanFret=cellfun(@(x,y) x/y,meanYFP,meanCFP,'UniformOutput',false);

% make DistMask
umPerPix=distanceScale(opt.objmag);
masks=cellfun(@(x) ~isnan(x)==1,imAcceptor,'UniformOutput',false);
masks=cellfun(@(x) returnMaskWithLargestArea(x),masks,'UniformOutput',false);
masks=cellfun(@(x) logical(x),masks,'UniformOutput',false);
% check if some mask images are empty. If empty flag this cell as bad and
% truncate all of the futher analysis.
cellMaskTF=cellfun(@(x) any(x,'all'),masks);
if cellMaskTF
    b.cellMaskPass=1;
     
    areaMask=cellfun(@(x) regionprops(x,'Area'),masks);
    distMasks=cellfun(@(x,y) (bwdistgeodesic(x,y,'quasi-euclidean')*umPerPix),masks,frapMask,'UniformOutput',false);
    
    %binMask=buildBinMasks(distMasks,opt.binedges,corrYFP,imDonor);
    
    binRange=determineBinRangesFromInputVect(opt.binedges);
    areaThresh=calcMinAreaThreshChunkyBin(distMasks,opt.numbefore,binRange);
    % need to end function if areaThresh is empty. areaThresh fails when the
    % frapspot is not in the cell.
    if isempty(areaThresh)
        b.areaThreshPass=0;
    else
        b.areaThreshPass=1;
        for i=1:length(distMasks)
            b(i).ratioCorYFP=corrYFP{i};
            b(i).binWidth=opt.binedges(2)-opt.binedges(1);
            b(i).distMask=distMasks{i};
            b(i).ratioBin0=wholeCellMeanFret{i};
            b(i).binArea0=areaMask(i).Area;
            b(i).areaThresh=areaThresh;
            for j=1:length(binRange)
                clear areaDat;
                tempBin=distMasks{i}>=binRange{j}(1) & distMasks{i}<binRange{j}(2);
                tempYFP=nanmean(corrYFP{i}(tempBin),'all');
                tempCFP=nanmean(imDonor{i}(tempBin),'all');
                tempRatio=tempYFP/tempCFP;
                
                b(i).(strcat('bin', num2str(j)))=tempBin;
                b(i).(strcat('YFPb', num2str(j)))=tempYFP;
                b(i).(strcat('CFPb', num2str(j)))=tempCFP;
                b(i).(strcat('ratioBin', num2str(j)))=tempRatio;
                %need to deal with empty ind especially for frame 1
                areaDat=regionprops(tempBin,'Area');
                if isempty(areaDat)
                    b(i).(strcat('binArea',num2str(j)))=nan;
                else
                    % some bins have multiple masks. I will sum the mask areas
                    sumArea=sum(arrayfun(@(x) x.Area,areaDat));
                    % use Bin 1 as min pix area thresh, store vals<thresh as nan
                    if sumArea >= areaThresh
                        b(i).(strcat('binArea',num2str(j)))=sumArea;
                    else
                        b(i).(strcat('binArea',num2str(j)))=nan;
                    end %if sumArea >=areaThresh
                end % if isempty(areaDat)
            end %for j length binRanges
        end % for i length  distMasks
    end %if empty area thersh.
    
else
    b.cellMaskPass=0;
    
end %if cellMaskTF
cBinOut=b;



end

% function binMask=buildBinMasks(distMask,binEdges,corrYFP,imDonor)
% binRange=determineBinRangesFromInputVect(binEdges);
% for i=1:length(distMask)
%     for j=1:length(binRange)
%         tempBin=distMask{i}>=binRange{j}(1) & distMask{i}<binRange{j}(2);
%         tempYFP=nanmean(corrYFP{i}(tempBin),'all');
%         tempCFP=nanmean(imDonor{i}(tempBin),'all');
%         tempRatio=tempYFP/tempCFP;
%         b(i).(strcat('bin', num2str(j)))=tempBin;
%         b(i).(strcat('YFPb', num2str(j)))=tempYFP;
%         b(i).(strcat('CFPb', num2str(j)))=tempCFP;
%         b(i).(strcat('ratioBin', num2str(j)))=tempRatio;
%         b(i).distMask=distMask{i};
%         
%     end
% end
% binMask=b;
% end

function bRange=determineBinRangesFromInputVect(binEdges)
bRange=cell(1,length(binEdges));
for r=1:length(bRange)
    if r==1
        bRange{r}=[0 binEdges(1)];
    else
        bRange{r}=[binEdges(r-1), binEdges(r)];
    end 
end
end