function GeneratePCAStalkSTL(stalknum,npts)
% FILENAME: GeneratePCAStalkSTL.m
% AUTHOR: Ryan Larson
% DATE: 1/31/2020
%
% PURPOSE: 
%   Generate an STL file that is an approximation of a specified stalk used
%   in LongPCAData.mat. Data process is loosely based on
%   GenerateStalkSegmentBetsy_V3.m.
% 
% 
% INPUTS:
%       
%       
% OUTPUTS:
%       - 
%
%
% NOTES: 
%       - 
% 
% 
% VERSION HISTORY:
% V1 - Assumes no meaningful drift (no xbar or ybar components included)
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load('LongPCAData.mat');

% Set up polar data structures
theta = zeros(length(keepcols),npts); % Rows are slices, columns are points around each cross-section
stalk_ext = zeros(length(keepcols),npts);
stalk_int = zeros(length(keepcols),npts);

% Verify that stalknum is part of keeprows. Throw an error if it's not a
% keeper.
if ~ismember(stalknum,keeprows)
    error('The chosen stalknum is not part of the longitudinal PCA set. Choose another stalk');
end

% Determine the adjusted index of the stalk of interest
stalkidx = find(keeprows == stalknum);

% Run PCA
[APCAs, Acoeffs, APCA_variances, Atstat, Aexplained, AvarMeans] = pca(APCA,'Centered',false);
[BPCAs, Bcoeffs, BPCA_variances, Btstat, Bexplained, BvarMeans] = pca(BPCA,'Centered',false);
[TPCAs, Tcoeffs, TPCA_variances, Ttstat, Texplained, TvarMeans] = pca(TPCA,'Centered',false);
[alphaPCAs, alphacoeffs, alphaPCA_variances, alphatstat, alphaexplained, alphavarMeans] = pca(alphaPCA,'Centered',false);

%% Populate polar data structures
% theta
for i = 1:size(theta,1)
    theta(i,:) = linspace(0,2*pi,npts);
end

% 


end
