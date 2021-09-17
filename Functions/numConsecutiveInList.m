function out=numConsecutiveInList(logicalVectIn)
%This function takes a logical vector as an input and counts the number of
%consecutive true values. The function will return a vector that is the
%same length of the input vector with >0 values indicating the number of
%consecutive logicals. The index for the >0 vals is the first position of
%the consecutive string.
out = double(diff([~logicalVectIn(1);logicalVectIn(:)]) == 1);
v = accumarray(cumsum(out).*logicalVectIn(:)+1,1);
out(out == 1) = v(2:end);
end