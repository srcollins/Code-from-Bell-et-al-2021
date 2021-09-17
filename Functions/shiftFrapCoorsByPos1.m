function adjFrapCoors = shiftFrapCoorsByPos1(frapCoors,pos1,cellIm,varargin)
%2020-07-16 GB: This function will return an nx2 matrix that contains the
%shifted frap coordinates for each image frame. The funtion uses the shift
%parameters in pos1 to make a list of new frap coors for each image in the
%experiment. The pos1 input should be a list of x,y coordinates that are in units
%of pixels. 

%2020-07-28 Update: the frap coordinates should be shifted only for the
%frames before stim when we are actually moving the stage.Additionally, all
%of the frap values should be in reference to frame 6.
%% process varargin inputs
opt.numbefore=5;


for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% duplicate the frapCoors into a matrix the length of 
adjFrapCoors=repmat(frapCoors,length(cellIm),1);
%% shift frap coors
for i=1:opt.numbefore
    
adjFrapCoors(i,:)=adjFrapCoors(i,:)+ [pos1(i,1)-pos1(6,1),pos1(i,2)-pos1(6,2)];

end