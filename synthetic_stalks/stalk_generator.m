% Ryan Larson - 3/14/2019
%
% Reformulation of GenerateStalkSegmentBetsy_v3.m using function methods,
% corrected equations for stalk cross section shape, and outputs data in an
% array that will be easy to work with for PCA of the longitudinal data.

clear
clc

%% Create data structure for holding lots of stalk information

% Use radius-z scheme for handling data so it doesn't get huge