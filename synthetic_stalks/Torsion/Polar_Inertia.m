function [polar_intertia] = Polar_Inertia(theta,rho,dr)
% FILENAME: Center_Inertia.m
% AUTHOR: Michael Ottesen
% DATE: 6/4/2020
%
% PURPOSE: Calculate the polar moment of inertia of a stalk cross section
% 
% INPUTS: 
%       theta: an array of points representing theta angle values
%
%       rho: an array of radius points for the desired shape
%
%       dr: Desired incriment step size moving out in the radius direction
%       
% OUTPUTS:       
%       polar_inertia: Total polar moment of inertia of the input shape
%
% NOTES: 
%       - 
% 
% PSEUDO-CODE:
%       Define number of points and segment size in theta and R direction
%       Begin nested loop to cylce each dr segment in each angle
%           Calculate radius to center of segment
%           Calculate polar moment of inertia for the segment
%       Sum up all the elements for total polar moment of inertia
%        
% VERSION HISTORY:
% V1 - Hard coded for shapes with constant dtheta and set dr = 0.1
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

% define number of points, dr (incriment size), and dtheta (angle incriment)
N = length(theta);

% Find the change in theta for each incriment in theta direction
for i = 1:N-1
    dtheta(i) = abs(abs(theta(i+1)) - abs(theta(i)));
end
dtheta(N) = abs(abs(theta(1))-abs(theta(N-1)));
% If there is a jump (ie. from almost 2*pi (like 6.1) to 0) for
% the last incriment, then add 2*pi to the low value to compensate.
if dtheta(N) > 3
    dtheta(N) = abs(abs(theta(1))+ 2*pi -abs(theta(N-1)));
else
end

% Cycle through all angles and dr value for each angle value 
for i = 1:N
    for j = 1:(length(0:dr:rho(i))-1)
        
        % Find radius length to center of segment, and polar moment of
        % inertia of small segment
        r = dr*j - (dr/2);
        J(i,j) = (r^3)*dtheta(i)*dr;
          
    end
end

% Find total polar inertia by summing all the small elements
polar_intertia = sum(J,'all');
end



