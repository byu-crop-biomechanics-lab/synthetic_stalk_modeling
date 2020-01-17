function AllTransversePCA(stalknums,slice_dists)
% FILENAME: AllTransversePCA.m
% AUTHOR: Ryan Larson
% DATE: 1/17/2020
%
% PURPOSE: 
% 
% 
% INPUTS:
%       
%       
% OUTPUTS:
%       
%
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
load StalksDCR_360pts.mat
hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

for slice = slice_dists
    slice
    % For each slice distance, iterate through the stalknums and get the
    % data to feed into the large PCA array
    
    % Gather all data for stalks at this slice distance
    
    % Flip cross-sections that need adjustment (only if there isn't flip
    % data). This covers all 980 stalks for the current slice distance
    
    % Get ellipse fits and difference data 
    
    % Add the difference data and corresponding ellipse fits to separate
    % large arrays for feeding into PCA. Take note of the starting and
    % ending indices corresponding to each slice distance so data can be
    % correctly reconstructed later. Also save the error_indices that
    % correspond.
    
    % Remove error cases and adjust indices, starting at the bottom of the
    % array. Save adjusted indices.
    
    
    % Run PCA on the resulting large data set. Save this for access by
    % transverse_wrapper_V4.m
    
    
    % MAKE SURE TO ADDRESS AND CHECK THE INTERPOLATION ISSUE THAT MIGHT
    % HAVE BEEN CAUSING ALL THE ERROR_INDICES AND THE NEED TO SHIFT
    % EVERYTHING
    
end


end