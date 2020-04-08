function [spline] = writespline_V2(len,data)
% FILENAME: writespline_V2.m
% AUTHOR: Ryan Larson
% DATE: 5/29/19
%
% PURPOSE: Turn vector spline data into a string
% 
% 
% INPUTS:
%       len: Length of the data (a holdover from previous versions, and
%       isn't fully necessary for good code)
% 
%       data: The original boundary vector
%       
% OUTPUTS:
%       spline: A string version of data that can be inserted in a Python
%       script
%
% NOTES:
%      
% 
% 
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
    for i = 1:len 
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), '); 
    end
end