function [maxind,minind] = boundinds(bound1,bound2,theta)
%boundinds.m: Find the indices in a vector theta that match bounds of
%interest

% Find indices for max and min theta values in range of interest:
inds = not(abs(sign(sign(bound1-theta) + sign(bound2-theta))));
maxind = max(find(inds,5,'last'));
minind = min(find(inds,5,'first'));

end

