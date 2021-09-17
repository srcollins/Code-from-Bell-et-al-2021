function drawShadedErrorRegionRMNans(x,y,errVals,varargin)
%2020-10-22 GB: If the input data is padded with nans the patch will fail.
%THis function will detect terminal nans and remove them.
%2021-03-18 GB: added facealpha as input.


%%process varargins
opt.color=[0.8 0.8 0.8];
opt.edgealpha=0.3;
opt.linestyle='--';
opt.facealpha=0.2;
for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% detect and remove nans
% check that x, y, and errVals are same length
if ~isempty(intersect3(numel(x), numel(y),numel(errVals)))
    xInd=~isnan(x);
    yInd=~isnan(y);
    evInd=~isnan(errVals);
else
    fprintf('numel of each input should be equal');
    return
end

% check that the indicies are equal

if isequal(xInd, yInd, evInd)
    
    x=x(xInd);
    y=y(yInd);
    errVals=errVals(evInd);
    
    patch([x(:); vect(x(end:-1:1))],[y(:)+errVals(:); vect(y(end:-1:1)-errVals(end:-1:1))],opt.color,...
        'FaceAlpha',opt.facealpha,'EdgeAlpha',opt.edgealpha,'HitTest','off','HandleVisibility','off','LineStyle', opt.linestyle);
else
    
    fprintf('numel of each input should be equal');
    return
end

end
