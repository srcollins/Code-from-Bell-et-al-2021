function ts=buildTSforPDLdataPreMetaFiles(folder,filenames1,varargin)
% 2020-10-14 GB: the purpose of this function is to create the time stamp
%(ts) data structure for PDL assay experiments before we were saving the 
%meta data files. I looked back at the experimental design for the pre meta
%data experiments and found that there are two protocol types. 

%type1: 1.5 sec interval, 10 frames before, numStim-1 between, 20 frames
%after. For stimulations stim1 was followed by pause(1). Then sequential
%imageing with TK then stim for the remainder. These images should be
%occurning within 1.5 sec. Assume that the time between TK and Dapi img is
%300ms.

%Type2: 1.5 sec interval, 10 frames before, 50 frames after. the stim is now part of the
%numafter for loop. Here the stim channel is imaged, pause(1), imgTK, pause
%for 1.5 sec interval remaining. THe good news is that we already
%calculated these times for the modeling. I will load timeData.mat and use
%these values to fill out ts.

%2020-11-11 GB Update to include detection of 2pulse data. THe processing
%method will change because the gaps between images are very different.


%2021-03-27 GB: AH and I noticed that time bars plotted with real meta data
%are not starting at t=0, instead its a few second delay. This function was
%developed using timing assumptions, however empirical data suggests that
%stimulation is slower. If the x axis is relative to the frame preceeding
%stimulation, then mean frap time is ~2seconds while the next tk frame is
%at 4 seconds. Update this function later.
%% process varargin
opt.stimterm={'DAPI','UV'};
opt.imgint=1.5;
opt.numbefore=10;
opt.datethresh='01-Sep-2019';

opt.int='(?<=Int_)[0-9]{0,3}';
opt.exp='(?<=Exp_)[0-9]{0,3}';
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% determine the number of stims
 %1) determine number of dapi or UV stims
    allImgFNames=getFilenames(folder,'img');
    numStims=[0 0];
    for i=1:length(opt.stimterm)
    numStims(i)=length(find(boolRegExp(allImgFNames,opt.stimterm{i})));
    end
    numStims=sum(numStims);
%% determine if exp was a 2 pulse protocol
%This will use the folder name as I was consistent in using the 2 pulse
%label
 twoPulse=boolRegExp(folder,'2Pulse');   
%% determine which stim protocol to use
% EXP type1=may 2018
%EXP type2=sept 2019
dateThresh=datetime(opt.datethresh);
imgDate=datestr(datetime(regexp(filenames1{1},'[\d]{8,8}(?=T)','match'),'Inputformat','yyyyMMdd'));
imgDate=datetime(imgDate);
if imgDate<dateThresh %Run TYPE1 TS:
    
    
    if twoPulse
        % NOTE: the 2 pulses are a bit of a pain because the 0 stims will
        % have timesgaps depending on the stimbetween parameter.
        secPostStim=str2double(regexp(folder,'(?<=secPostStim1_).*','match')); % these vals are in rounded seconds. Its easier to conver to frame number
        framesBetweenStim=round(secPostStim/opt.imgint);
        % get stimPwer
        int=regexp(folder,opt.int,'match');
        exp=regexp(folder,opt.exp,'match');
        pwr=str2double(int)*str2double(exp);
        tkTime=[0:opt.imgint:(length(filenames1)-1)*opt.imgint];
        
        if pwr==0
            % for no stim there was a 1 sec pause, but no delay for image
            % capture. add stim 2 pause use 1 sec.
            stimPause=1;
            stim2Ind=opt.numbefore+framesBetweenStim+1; % need to add 1 because framesBetween is the last TK frame Before second Stim
            tkTime(stim2Ind:end)=tkTime(stim2Ind:end)+stimPause;
            % add stim1
            tkTime(opt.numbefore+1:end)=tkTime(opt.numbefore+1:end)+stimPause;
            stimTime=nan;
        else
            %add stim 2 pause use 1.3 sec. 1 sec is the pause, 0.3 for img snap
            stimPause=1.3;
            stim2Ind=opt.numbefore+framesBetweenStim+1; % need to add 1 because framesBetween is the last TK frame Before second Stim
            tkTime(stim2Ind:end)=tkTime(stim2Ind:end)+stimPause;
            % add stim1
            tkTime(opt.numbefore+1:end)=tkTime(opt.numbefore+1:end)+stimPause;
            stimTime=[0 0];
            stimTime(1)=tkTime(opt.numbefore)+.3;
            stimTime(2)=tkTime(opt.numbefore+framesBetweenStim)+0.3;
        end
        stimTime=stimTime-tkTime(opt.numbefore);
        tkTime=tkTime-tkTime(opt.numbefore);
        
       
        
    else
        if numStims==0
            tkTime=[0:opt.imgint:(length(filenames1)-1)*opt.imgint];
            tkTime=tkTime-tkTime(opt.numbefore);
            stimTime=nan;
        elseif numStims==1
            tkTime=[0:opt.imgint:(length(filenames1)-1)*opt.imgint]; % get time vect for all img
            tkTime(opt.numbefore+1:end)=tkTime(opt.numbefore+1:end)+1.3; % add in the 1.3 sec pause for stim 1
            tkTime=tkTime-tkTime(opt.numbefore); % subtract time at numbefore to set the time for frame 10 at 0
            stimTime=tkTime(opt.numbefore)+0.3; % Define the stim time relative to the tkTim
        elseif numStims>1 
            tkTime=[0:opt.imgint:(length(filenames1)-1)*opt.imgint]; % get time vect for all img
            tkTime(opt.numbefore+1:end)=tkTime(opt.numbefore+1:end)+1.3; % add in the 1.3 sec pause for stim 1
            tkTime=tkTime-tkTime(opt.numbefore);
            stimTime=zeros(1,numStims);
            for i=1:length(stimTime)
                fRange=[opt.numbefore:opt.numbefore+length(stimTime)];
                stimTime(i)=tkTime(fRange(i))+0.3;
            end
            
        end
    end %if twoPulse
    
else
    %Run Type 2 TS
    load([folder filesep 'timeData.mat']);
    tkTime=imTimePost;
    tkTime=tkTime-imTimePost(opt.numbefore);
    stimTime=stimTimePost;
    stimTime=stimTime-imTimePost(opt.numbefore);
    
end  %if imgDate<dateThresh

ts.tkTime=tkTime;
ts.frapTime=stimTime;
ts.frameRate=round(mean(diff(ts.tkTime(1:10))),2);
ts.twoPulse=twoPulse;

end