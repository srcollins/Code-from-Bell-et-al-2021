function [yOut,numOut] = evalBySlidingBins(xIn,yIn,bins,binwidth,func)
%syntax: res = evalBySlidingBins(xIn,yIn,bins,func)  Evaluates the function
%specified by func on the data (xIn,yIn) by binning the x values. It bins
%the x values using sliding bins centered at the locations specified by
%bins and with width binwidth. The function func should take a list and
%return a value.

x2=xIn(:); y2=yIn(:);
tempInd= isnan(x2) | isnan(y2);
x3=x2(~tempInd);
y3=y2(~tempInd);
mat=[x3 y3];
smat=sortrows(mat,1);
b1=1; b2=1;
yOut=nan(size(bins));
numOut=nan(size(bins));
for i=1:length(bins)
%     b1=findmin(smat(:,1),bins(i),binwidth,b1);
%     b2=findmax(smat(:,1),bins(i),binwidth,b2);
    b1=find(smat(:,1)>=bins(i)-binwidth,1,'first');
    b2=find(smat(:,1)<bins(i)+binwidth,1,'last');
    list=smat(b1:b2,2);
    if ~isempty(list)
        yOut(i)=feval(func,list);
        numOut(i)=length(list);
    end
end

end

%*************************************************************************

function ind=findmin(list,middle,binsize,prev)
target=middle-binsize/2;
ind=prev;
while (list(ind)<target)&(ind<length(list))
    ind=ind+1;
end
end

%*************************************************************************

function ind=findmax(list,middle,binsize,prev)
target=middle+binsize/2;
ind=prev;
while (list(ind)<target)&(ind<length(list))
    ind=ind+1;
end
end

