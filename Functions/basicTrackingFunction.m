function traj=basicTrackingFunction(coord,maxDisp)
% A streamlined tracking function that uses annsearch to identify nearest
% neighbors, and then selects reciprocal best matches between adjacent
% frames.
% coord should be a cell array of cell coordinates
% maxDisp should be the maximum allowed displacement in pixels

trajIndMat=nan(sum(cellfun(@(x) size(x,1),coord(1:end-1))),length(coord));
trajIndMat(1:size(coord{1},1),1)=vect(1:size(coord{1},1));

numTraj=size(coord{1},1);
for i=1:(length(coord)-1)
    if ~isempty(coord{i}) && ~isempty(coord{i+1})
        % Forward match
        [ix,d]=annsearch(coord{i+1}(:,1:2)',coord{i}(:,1:2)',1);
        M1=[vect(1:length(ix)) ix(:)];
        M1(d>maxDisp,:)=[];
        
        % Backwards match
        [ix,d]=annsearch(coord{i}(:,1:2)',coord{i+1}(:,1:2)',1);
        M2=[vect(1:length(ix)) ix(:)];
        M2(d>maxDisp,:)=[];
        
        % Filter for reciprocal matches
        M=intersect(M1(:,1:2),M2(:,[2 1]),'rows');
    else
        M=[0 0];   % There can't be any matches if there are no cells in one of the frames. [0 0] is used as a placeholder to avoid an error from intersect
    end
    
    % Apply matches
    [junk,i1,i2]=intersect(trajIndMat(:,i),M(:,1));
    trajIndMat(i1,i+1)=M(i2,2);   %Continuing trajectories
    if i<(length(coord)-1)      %Only need to add new trajectories if it is not the last frame
        idVect=vect(1:size(coord{i+1},1));
        [junk,i1]=setdiff(idVect,M(:,2));
        trajIndMat((numTraj+1):(numTraj+length(i1)),i+1)=idVect(i1);  %Potential new trajectories
        numTraj=numTraj+length(i1);
    end
end

% Clean up trajIndMat
trajIndMat=trajIndMat(1:numTraj,:);   %Remove empty rows
hasAtLeastOneLink=sum(~isnan(trajIndMat),2)>1;
trajIndMat=trajIndMat(hasAtLeastOneLink,:);

% Generate a cell array of trajectories
numTraj=size(trajIndMat,1);
numCols=size(coord{1},2)+3;
traj=cell(numTraj,1);
for i=1:numTraj
    framesCovered=find(~isnan(trajIndMat(i,:)));
    traj{i}=nan(length(framesCovered),numCols);
    for j=1:length(framesCovered)
        frameNum=framesCovered(j);
        traj{i}(j,:)=[coord{frameNum}(trajIndMat(i,frameNum),:) trajIndMat(i,frameNum) i frameNum];
    end
end
