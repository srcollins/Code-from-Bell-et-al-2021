function ax=density_scatter_heatmap(v1,v2,nbins1,nbins2,ax,scale)

if nargin<3
    nbins1=50;
    nbins2=50;
end

% Make sure that v1 and v2 are vectors
v1=v1(:);
v2=v2(:);

% Range for v1
if length(nbins1)==1
    min1=min(v1);
    max1=max(v1);
    step1=(max1 - min1)/nbins1;
    edges1=min1:step1:max1;
    ctrs{1}=0.5*(edges1(1:(end-1)) + edges1(2:end));
else
    ctrs{1}=nbins1;     % nbins1 was entered as an array of centers rather than as a number of bins
end

% Range for v2
if length(nbins2)==1
    min2=min(v2);
    max2=max(v2);
    step2=(max2 - min2)/nbins2;
    edges2=min2:step2:max2;
    ctrs{2}=0.5*(edges2(1:(end-1)) + edges2(2:end));
else
    ctrs{2}=nbins2;
end

mat=hist3([v1 v2],ctrs);
if nargin<6 | strcmpi(scale,'log')
    mat=log(mat+0.5);   % The color will be determined by the log density of points
else
    %linear scale
end
if nargin<5 || isempty(ax)
    figure;
else
    axes(ax);
end
imagesc(ctrs{1},ctrs{2},mat');

ax=gca;
set(ax,'YDir','normal');
set(ax,'box','off');
colormap(parula_black);
