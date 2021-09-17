function bleachedMol=bleachedMoleculeDiffusionFunction(p,x,t)
% This function returns the concentration profile of bleached sensor
% molecules after a local laser pulse at time zero. It assumes a diffusion
% coefficient of 0.5 microns squared per second, and an empirically fit
% initial bleaching pattern
%
% Inputs:
% p = parameters. p(1) is a magnitude parameter which should be
% proportional to laser output power. p(2) is a shape parameter which
% describes the initial bleaching pattern.
% x is the distances from the center of the laser spot
% t is time since the laser pulse

% The function is made from an analytic solution to the radial diffusion
% problem.
if p(1)==0 || t<0
    bleachedMol=ones(size(x));
else
    D=0.5; % diffusion coefficient in microns squared per second
    bleachedMol = p(1)*exp((-x.^2)/(D*4*(t+p(2)))) / (4*pi*D*(t+p(2)));
end


end