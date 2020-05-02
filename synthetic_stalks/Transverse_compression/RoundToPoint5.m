function [roundval] = RoundToPoint5(input)
% FILENAME: getBoxLims.m
% AUTHOR: Ryan Larson
% DATE: 4/22/2020
%
% PURPOSE: Round a number to a the closest 0.5
% 
% INPUTS:
%       input: A single float value
%       
% OUTPUTS: 
%       roundval: The output value, rounded to the nearest 0.5.
%
% NOTES: 
% 
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Calculate the remainder when dividing the input number by 0.5.
%   If the remainder is greater than or equal to 0.25:
%       Round up to the next 0.5
%   else:
%       Round down to the next 0.5
%   end
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

% Calculate the remainder when dividing the input number by 0.5
modval = mod(input,0.5);

% If the remainder is greater than or equal to 0.25, round up to the next
% 0.5. Otherwise, round down to the next 0.5.
if modval >= 0.25
    roundval = ceil(input*2)/2;
else
    roundval = floor(input*2)/2;
end

end