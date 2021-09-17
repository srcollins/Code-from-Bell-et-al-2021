function frapInd = convertXYtoRC(ind)
%this function will convert frap coordinates from XY position to Row by
%column position or vice versa. In matlab, matricies are organized in a row by column
%format, where rows are equivalent to the Y axis and the X axis is
%equivalent to the number of columns.

%Input:[X,Y];
%%
frapInd=nan(1,2);
reqDim=[1 2];
if size(ind)==reqDim
frapInd(1)=ind(2);
frapInd(2)=ind(1);
else
    fprintf('input must be a 1 by 2 vector');
    return
end
