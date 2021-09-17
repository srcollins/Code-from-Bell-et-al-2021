function [p0] = buildP0TomKat(varargin)
%% This function builds p0, an input for computeRegistrationParametersQuadFit
%that contains the user's best guess for the initial translation alignment
%parameter. The function when called with no inputs will default to good
%TomKat alignment guess. Use the name value pairs to add in your own
%inputs.

%% read in the varargins
opt.x1=-10;
opt.y1=-39;

for i=1:2:length(varargin)
    opt.(lower(varargin{i}))=varargin{i+1};
end
%% build p0
 p0.pX=zeros(6,1); p0.pX(1)=opt.x1;
 p0.pY=zeros(6,1); p0.pY(1)=opt.y1;
end

