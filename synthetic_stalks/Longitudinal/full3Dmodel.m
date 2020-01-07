function full3Dmodel(DataTable,totstalks,nstalks)

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
%       DataTable - Table with A, B, T, and theta_rot values, as well as
%       cross-section longitudinal position and stalk number (simplified
%       from the original downsampled data table, like the DCR table was
%       
%       totstalks - Total number of stalks used for PCA (should be close to
%       980)
%       
%       nstalks - Number of stalks to create 
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

% Get rid of problem cases where ellipse fitting doesn't work

% Make sure notches are all pointed the correct way (might need a
% notch-direction checking function to do this)


%% Select stalks

% Read the input number of stalks and get a random, ordered array of stalk
% numbers
stalks = randperm(980,nstalks);
stalks = sort(stalks);



%% Create PCA data table if not already created
% Check if PCA variable table exists in the current directory. If yes, skip
% this section

PCAvars = NaN(totstalks,(nslices*4)); % 4 variables: A,B,T,theta_rot

% Fill in PCAvars with the corresponding A,B,T,theta_rot data

% Run PCA on PCAvars table
% Save output table as .mat file


%% Create geometric cases and corresponding Python scripts


end