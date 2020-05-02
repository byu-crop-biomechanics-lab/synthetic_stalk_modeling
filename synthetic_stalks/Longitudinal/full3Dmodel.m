function full3Dmodel(stalknums,DataTable,totstalks,nstalks)

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
%       stalknums - A vector of unique integers from 1 to 980 that determines
%       which stalks to sample from (use randperm(980,K) to choose K
%       unique integers from 1 to 980). Sort this vector so stalk numbers
%       are in ascending order
% 
%       DataTable - Table with A, B, T, and theta_rot values, as well as
%       cross-section longitudinal position and stalk number (simplified
%       from the original downsampled data table, like the DCR table for
%       transverse deformation
%       
%       totstalks - Total number of stalks used for PCA (should be close to
%       980) MIGHT NOT NEED THIS INPUT
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