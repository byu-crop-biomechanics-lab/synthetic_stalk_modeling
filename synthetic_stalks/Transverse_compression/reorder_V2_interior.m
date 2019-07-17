function [X, Y, x, y, r, t, xp, yp, rp, tp] = reorder_V2_interior(x, y, alpha, ordx, ordy)
% function centers and reorders the parametric curves described
% by x and y.
%
% INPUTS:   x,y - column vectors of x & y values
%           alpha - angular orientation of the stalk  
%           ordx,ordy - xy values for the center of rotation of the
%           exterior
%
% OUTPUTS:  Nomenclature: CAPS variables in image coordinates, lower case relative to stalk coordinate system. "p" indicates "prime" - for a rotated coordinate system.
%           X and Y - reordered x and y values.
%           x, y - reordered x and y values relative to xbar and ybar (but not rotated)
%           r, t - reordered polar coordinates: r and theta values (NOT rotated)
%           xp, yp - reordered x and y values relative to xbar and ybar AND rotated.
%           rp, tp - reordered polar coordinates: r and theta values, AND rotated
%
% Note: designed for use with external contour only (hard to do both since int and ext have different numbers of data points).
%
% VERSION HISTORY
%   V2 - 4/5/18 - updated inputs and comments to reflect changes since conversion from smorder_R2
%
%   V4 - 7/17/2019 - Made version with no shift inside the function, so
%   shifting can take place externally. Allows interior and exterior
%   boundaries to be shifted identically outside the reorder_V4 function.

% 1. shift origin to center
X = x;
Y = y;

% Locate Center
xbar = ordx;
ybar = ordy;

% shift origin to center
y = y - ybar;
x = x - xbar;

% 2. compute theta values in the (xbar, ybar) coordinate system
t = atan2(y,x);
A = size(t);
if A(1) == 1
    t = t';
    x = x';
    y = y';
    X = X';
    Y = Y';
end

% 3. correct progression direction    
% make sure points progress in the positive direction - this is done by checking
% the sign of the median difference between theta values. If positive, rotation is positive. If
% negative, all indices are reversed.
theta_diff = diff(t);                   % differences between subsequent theta values
med_diff = median(theta_diff);              % find the MEDIAN theta difference (required because a theta shift value of +/-2pi is always present)

if med_diff < 0                             % check if theta differences are negative
    t = t(end:-1:1);                % if so, reverse indices in all variables
    x = x(end:-1:1);
    y = y(end:-1:1);
    X = X(end:-1:1);
    Y = Y(end:-1:1);
end

% 4. Reorder variables starting from the value closest to alpha and shift theta values to the stalk coordinate system
temp_theta = t;                         % temp_theta variable used for convenience
index = (t<=alpha);                     % indices of theta values less than alpha
temp_theta(index) = NaN;                    % these values are ignored using NaN
[minval, alphloc] = min(temp_theta);         % find the smallest remaining value
% size(t(alphloc:end))
% size(t(1:alphloc-1))
t = [t(alphloc:end); t(1:alphloc-1)];      % all indices are reordered accordingly.
x = [x(alphloc:end); x(1:alphloc-1)];
y = [y(alphloc:end); y(1:alphloc-1)];
X = [X(alphloc:end); X(1:alphloc-1)];
Y = [Y(alphloc:end); Y(1:alphloc-1)];

% x = x + xbar;
% y = y + ybar;

% 5. Shift theta values so that they are continuous (no jumps) 
% this loop is fairly wasteful, costs about 0.04 seconds.  There are much faster ways, for example just searching for sudden jumps in theta.  But these are likely less robust.
for i = 2:length(t)    
    poss_theta = [t(i) - 2*pi, t(i)-pi, t(i), t(i) + pi, t(i)+2*pi];
    [minval, alphloc] = min(abs(t(i-1) - poss_theta));
    t(i) = poss_theta(alphloc);   
end
% IS THIS STEP EVEN NECESSARY????


% 6. Define radius values
rsq = x.^2 + y.^2;


% 7. Data for coordinate system centered at (xbar, ybar) (no tilt) - First point is at the major axis
x = x; 
y = y;
r = rsq.^0.5;
t = t;

% 8. Data for tilted coordinate system centered at (xbar, ybar) (no tilt)
xp = x*cos(-alpha)-y*sin(-alpha);
yp = x*sin(-alpha)+y*cos(-alpha);
rp = (x.^2 + y.^2).^0.5;
tp = t-alpha;




