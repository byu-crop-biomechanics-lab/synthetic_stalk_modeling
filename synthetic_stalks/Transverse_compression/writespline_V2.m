function [spline] = writespline_V2(data)
% FILENAME: writespline_V2.m
% AUTHOR: Ryan Larson
% DATE: 5/29/19
%
% PURPOSE: Turn spline data from a 2-column array into a string
% 
% 
% INPUTS:
%       data: The original boundary array (2 columns, where column 1 is x
%       data and column 2 is y data)
%       
% OUTPUTS:
%       spline: A string version of data that can be inserted in a Python
%       script
%
% NOTES:
%      
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Create an empty string to hold the spline data (or the ordered XY pairs
%   that define the exterior or interior cross-section profile).
% 
%   For each row in the data (a 2-column array of XY data):
%       Concatenate the current data point to the spline in the ordered
%       pair form, like (x,y).
%   end
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - Writes spline to a text file which can then be copied manually into 
% V2 - Made writespline a function that works with existing functions
% instead of writing the spline to a text file
% V3 - 
%
% -------------------------------------------------------------------------
    %define empty spline and number of x-y points
    spline = '';
    
    %run through 1-column arrays of the x and y data points for the spline, and add to the end of the string with the correct formatting
    for i = 1:length(data)
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), ');
    end
end