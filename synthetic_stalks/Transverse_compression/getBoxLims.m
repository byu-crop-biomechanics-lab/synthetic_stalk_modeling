function [uplim,lolim] = getBoxLims(datacol)
% FILENAME: getBoxLims.m
% AUTHOR: Ryan Larson
% DATE: 4/22/2020
%
% PURPOSE: Take in a column vector of data going into a box plot and 
% determine reasonable upper and lower limits for the whiskers. These won't
% be exactly what Matlab puts out, but it will be close.
% 
% 
% INPUTS:
%       datacol: A column vector with the data that will go into the box
%       plot
%       
% OUTPUTS: 
%       uplim: The expected upper limit for the box plot whiskers.
% 
%       lolim: The expected lower limit for the box plot whiskers.
%
% NOTES: 
% 
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Calculate the first and third quantiles.
%   Calculate the interquartile range.
%   Calculate the upper and lower whisker positions by adding or
%   subtracting 1.5 times the interquartile range 
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - Simple estimation of whisker locations in a box plot
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

Q1 = quantile(datacol,0.25);
Q3 = quantile(datacol,0.75);
IQR = Q3 - Q1;

uplim = Q3 + 1.5*IQR;
lolim = Q1 - 1.5*IQR;

end