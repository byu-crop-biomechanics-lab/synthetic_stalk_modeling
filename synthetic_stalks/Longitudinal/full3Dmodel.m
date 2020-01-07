% FILENAME: full3Dmodel.m
% AUTHOR: Ryan Larson
% DATE: 1/7/2020
%
% PURPOSE: Take the methods from transverse_wrapper_V3.m, fix the problems
% with indexing, and apply the basic methods to longitudinal PCA to create
% a full 3D model of a stalk
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PSEUDO-CODE %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read from table containing original data (This might need to be a separate prepping function)

% Do stuff from table preparing function PrepSections_V2.m, except store
% the rotations of the cross-sections for later (keep original detected
% values, as well as relative values with respect to the rotation at the
% node

% Downsample each cross-section so there's less data to deal with? Maybe
% this isn't necessary?

% Rotate each cross-section in a given stalk so the rotation of the major
% axis of the node cross-section is 0 degrees from the x-axis and the notch
% is on the left for cross-sections below the node


%% 





