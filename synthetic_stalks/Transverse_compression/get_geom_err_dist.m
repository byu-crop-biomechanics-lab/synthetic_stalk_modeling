function [geom_err_dist] = get_geom_err_dist(ChosenEllipseData,PCAData,NewGoodStalks,nNEPCs)
% Correct issues in geometric_error.m so a set of stalks can be examined
% for geometric accuracy in a distribution format, and the results combined
% across longitudinal sampling location.

hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

load(PCAData,'ext_rhoPCAs','ext_rhocoeffs');
load(ChosenEllipseData,'ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int','R_ext','R_int','AVG_RIND_T','B');
load(NewGoodStalks,'newgoodstalknums');

% Create output data structure
M = length(AVG_RIND_T)*nNEPCs;
N = nNEPCs + 1;
geom_err_dist = zeros(M,N);


%% Loop
% For each stalk, take the difference between the original exterior data
% and the current geometric case (ellipse + some NEPCs).

for i = 1:M
    real_ext = R_ext(i,:);
    ellipse_ext = ELLIPSE_R_ext(i,:);
    
    % Calculate the differences for the ellipse-only case
    geom_err = zeros(size(real_ext));
    for j = 1:length(geom_err)
        geom_err(j) = abs(ellipse_ext(j) - real_ext(j)); % Should this be signed error or magnitude only?
    end
    
    % Take the sum of the error around the cross-section and insert it into
    % the output structure in the first column
    tot_err = sum(geom_err);
    geom_err_dist(i,1) = tot_err;
    
    for j = 1:N-1
        % Calculate the geometry of the current NEPC case
        [NEPC_ext,~] = section_from_PCA(newgoodstalknums,i,ELLIPSE_R_ext,ext_rhoPCAs,ext_rhocoeffs,AVG_RIND_T,j);
%         polarplot(ELLIPSE_T(i,:),NEPC_ext);
%         pause();
        
        % Calculate the differences for the current NEPC case
        geom_err = zeros(size(real_ext));
        for k = 1:length(geom_err)
            geom_err(k) = abs(NEPC_ext(k) - real_ext(k)); % Should this be signed error or magnitude only?
        end
        
        % Take the sum of the error around the cross-section and insert it into
        % the output structure in the first column
        tot_err = sum(geom_err);
        geom_err_dist(i,j+1) = tot_err;

    end
    
end

    

end