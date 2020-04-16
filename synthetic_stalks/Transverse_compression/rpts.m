function [r] = rpts(N,theta,dmaj,dmin)
% FILENAME: rpts.m
% AUTHOR: Ryan Larson
% DATE: 12/4/19
%
% PURPOSE: Simple function to turn ellipse diameter parameters into polar
% values, given a desired number of points.
% 
% 
% INPUTS:
%       N: The number of points in theta (redundant but haven't gone back
%       to change this in everything downstream)
% 
%       theta: Linearly-spaced vector of theta values from 0 to 2*pi.
% 
%       dmaj: Major diameter of the ellipse.
% 
%       dmin: Minor diameter of the ellipse.
%       
% OUTPUTS:
%       r: The radial value of the ellipse at each value of theta.
%
% NOTES:
%       
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Create an empty row vector of radius values with length N.
% 
%   for each element of the empty row vector:
%       Calculate the radial value of the ellipse, given the major and
%       minor diameters and the current angle.
%   end
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end