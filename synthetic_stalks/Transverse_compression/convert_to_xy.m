function [xy_columns] = convert_to_xy(R,theta)
% FILENAME: writespline_V2.m
% AUTHOR: Ryan Larson
% DATE: 5/29/19
%
% PURPOSE: Converts from polar to Cartesian
% 
% 
% INPUTS:
%       R: Radius data vector
% 
%       theta: Angle vector
%       
% OUTPUTS:
%       xy_columns: A 2-column array of xy data describing the previously
%       polar data
%
% NOTES:
%      
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
    N = length(theta);
    xy_columns = zeros(N,2);
    for i = 1:N
        xy_columns(i,1) = R(i)*cos(theta(i));
        xy_columns(i,2) = R(i)*sin(theta(i));
    end
end