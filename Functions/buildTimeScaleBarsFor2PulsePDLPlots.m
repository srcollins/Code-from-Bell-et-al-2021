function timeBars = buildTimeScaleBarsFor2PulsePDLPlots(pDatPDL,minY,inds2plot, colors, varargin)
%2020-10-21 GB: This function is designed to add stimulation durtation
%scalebars to PDL assay plots. These bars should also match the color of
%the corresponding curve.

%2020-12-15 GB: This function was updated and renamed to buildTimeScaleBarsFor2PulsePDLPlots
%so that it would automatically plot gray scale bars for the two pulse PDL
%plots

%2021-03-28 GB Update: added varargin inputs to use fieldnames in pDatPDL
%other than meanStimTime and meanTKtime; 

%% process varargin 


%%
opt.yaxismax=1.07;
opt.stimfieldname='meanStimTime';
opt.tktimefieldname='meanTKtime';
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% make time bars
% determine if any of the inds2plot are 0 stim conditions.
stimInd=[];
for i=1:length(inds2plot)
    stimInd(i)=any(~isnan(pDatPDL(inds2plot(i)).(opt.stimfieldname)));
end

%preallocate timebars
for j=sum(stimInd):-1:1
   timeBars(j).x=[];
   timeBars(j).y=[];
   timeBars(j).w=[];
   timeBars(j).h=[];
   timeBars(j).colors=[];
end

%% determine placement
inds2PFinal=inds2plot(logical(stimInd));

%apply stimInds to colors
stimIndColors=logical(stimInd');
%stimIndColors=logical(repmat(stimIndColors,1,3));
colors(~stimIndColors,:)=[];

% fill in timeBars
count=0;
for t=1:length(inds2PFinal)
    nextTKimgInd=[];
    nextTKimg=[];
    timeBars(t).x=pDatPDL(inds2PFinal(t)).(opt.stimfieldname);
    % calculate width for each pulse
    for j=1:length(timeBars(t).x)

        nextTKimgInd=find(pDatPDL(inds2PFinal(t)).(opt.tktimefieldname) > timeBars(t).x(j));
        nextTKimg=pDatPDL(inds2PFinal(t)).(opt.tktimefieldname)(nextTKimgInd(1));
        timeBars(t).w(j)=nextTKimg-timeBars(t).x(j);
        timeBars(t).y(j)=minY;
        timeBars(t).h(j)=opt.yaxismax-minY;
    end
    timeBars(t).colors=colors(t,:);
end

end