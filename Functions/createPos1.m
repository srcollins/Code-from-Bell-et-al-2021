function [pos1, pos1PassFlag]=createPos1(folder,objMag,rotationMatrix, varargin)
%% 2020-07-16 GB: This function creates Pos1 a list of pixel shifts based
%on stage movement

% process varargin
opt.filename='stage_positions.txt';

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end

%% create pos1
fid=fopen([folder filesep opt.filename],'r');
lines=readAllLines(fid);
if length(lines)>1
    fclose(fid);
    pos1=nan(length(lines),2);
    pos1(:,1)=str2double(getTokens(lines,'x\=([\-0-9\.]+)\,'));
    pos1(:,1)=pos1(:,1)-pos1(2,1);
    pos1(:,2)=-1*str2double(getTokens(lines,'y\=([\-0-9\.]+)'));
    pos1(:,2)=pos1(:,2)-pos1(2,2);
    pos0=pos1;
    % Convert stage shifts in microns to image shifts in pixels
    pos1=pos0/distanceScale(objMag);
    pos1=(rotationMatrix*pos1')';
    pos1=round(pos1);
    %remove nans from pos1 that are introduced by the multiStim lines
    indPos=isnan(pos1(:,1));
    pos1(indPos,:)=[];
    pos1PassFlag=true;
else
    pos1PassFlag=false;
    pos1=nan(1,2);
end
end