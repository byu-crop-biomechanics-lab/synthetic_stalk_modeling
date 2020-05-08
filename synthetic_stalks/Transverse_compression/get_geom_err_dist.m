function [geom_err_dist] = get_geom_err_dist(ChosenEllipseData,PCAData,NewGoodStalks,nNEPCs)
% FILENAME: get_geom_err_dist.m
% AUTHOR: Ryan Larson
% DATE: 2019
%
% PURPOSE: Gather geometric error data for each cross-section under
% examination. The inputs are somewhat limiting in that they only allow
% examination of a single slice location at a time (i.e. 5mm above the node
% only).
% 
% 
% INPUTS:
%       ChosenEllipseData: .mat file with the ellipse fit data for the
%       chosen slices.
% 
%       PCAData: Principal component .mat file for the chosen slices.
% 
%       NewGoodStalks: .mat file containing the the stalk numbers that were
%       used for the data set (? This header was written long after I had 
%       moved on to another version of the code).
% 
%       nNEPCs: The number of principal components to include in the
%       output.
%       
% OUTPUTS:
%       geom_err_dist: An array where rows correspond to unique
%       cross-section slices, and the columns correspond to different
%       levels of geometric approximation. First column is the ellipse
%       case, and the subsequent columns are the principal components in
%       numerical order. The content of the array is the radial error for
%       each theta coordinate, for each unique cross-section in the data.
% 
% NOTES: 
%       Follow up with plot_intervals_geom_err.m to visualize the
%       progression of geometric error as more principal components are
%       included in the model
% 
% -------------------------------------------------------------------------
% SUBROUTINES:
%   
% 
% PSEUDO-CODE:
%   
%       
% -------------------------------------------------------------------------
% 
% VERSION HISTORY:
% V1 - Correcting issues from geometric_error.m so a set of stalks can be
% examined for geometric accuracy in a distribution format, and the results
% manually combined across longitudinal sampling locations.
% V2 - 
% V3 - 
% -------------------------------------------------------------------------

hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

load(PCAData,'ext_rhoPCAs','ext_rhocoeffs');
load(ChosenEllipseData,'ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int','R_ext','R_int','AVG_RIND_T','B');
load(NewGoodStalks,'newgoodstalknums');

% Create output data structure
npts = 360;
M = length(newgoodstalknums)*npts;
m = length(newgoodstalknums);
N = nNEPCs + 1;
geom_err_dist = zeros(M,N);


%% Loop
% For each stalk, take the difference between the original exterior data
% and the current geometric case (ellipse + some NEPCs).

for i = 1:m
    starting_ind = (i-1)*npts + 1;    
    
    real_ext = R_ext(i,:);
    ellipse_ext = ELLIPSE_R_ext(i,:);
    
    % Calculate the differences for the ellipse-only case
    geom_err = zeros(size(real_ext'));
    for j = 1:length(geom_err)
        geom_err(j) = 100*(ellipse_ext(j) - real_ext(j))/B(i);
    end
    
    % Place the ellipse errors in the array
    geom_err_dist(starting_ind:(starting_ind+npts-1),1) = geom_err;
    
    for j = 1:N-1
        % Calculate the geometry of the current NEPC case (i is the index
        % in newgoodstalknums)
        [NEPC_ext,~] = section_from_PCA(newgoodstalknums,i,ELLIPSE_R_ext,ext_rhoPCAs,ext_rhocoeffs,AVG_RIND_T,j);
%         polarplot(ELLIPSE_T(i,:),NEPC_ext);
%         pause();
        
        % Calculate the differences for the current NEPC case
        geom_err = zeros(size(real_ext'));
        for k = 1:length(geom_err)
            geom_err(k) = 100*(NEPC_ext(k) - real_ext(k))/B(i); % Should this be signed error or magnitude only?
        end
        
        % Insert the calculated differences in the error array
        geom_err_dist(starting_ind:(starting_ind+npts-1),j+1) = geom_err;

    end
    
end

    

end