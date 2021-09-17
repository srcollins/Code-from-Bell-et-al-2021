function xMat=getXMatGridForImageRegistrationQuadFit(xvals,yvals)

[x,y]=meshgrid(xvals,yvals);
x=x(:); y=y(:);
xMat=ones(length(x(:)),6);
xMat(:,2)=x;
xMat(:,3)=y;
xMat(:,4)=x.*x;
xMat(:,5)=y.*y;
xMat(:,6)=x.*y;

