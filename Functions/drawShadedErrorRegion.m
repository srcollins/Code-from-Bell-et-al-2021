function drawShadedErrorRegionRMNans(x,y,errVals,c)

if nargin<4
    c=[0.8 0.8 0.8];
end
patch([x(:); vect(x(end:-1:1))],[y(:)+errVals(:); vect(y(end:-1:1)-errVals(end:-1:1))],c,'FaceAlpha',0.2,'EdgeAlpha',0,'HitTest','off','HandleVisibility','off');