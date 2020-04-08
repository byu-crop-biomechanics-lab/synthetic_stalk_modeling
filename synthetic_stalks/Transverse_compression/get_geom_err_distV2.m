function [geom_err_dist_ext,geom_err_dist_int] = get_geom_err_distV2(AllSlicesPCA,numPCs)
% Correct issues in geometric_error.m so a set of stalks can be examined
% for geometric accuracy in a distribution format, and the results combined
% across longitudinal sampling location.

hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

load(AllSlicesPCA);

% load(PCAData,'ext_rhoPCAs','ext_rhocoeffs');
% load(ChosenEllipseData,'ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int','R_ext','R_int','AVG_RIND_T','B');
% load(NewGoodStalks,'newgoodstalknums');

% Create output data structure
npts = 360;
M = size(ALL_DIFF_R_ext,1)*npts;
m = size(ALL_DIFF_R_ext,1);
N = numPCs + 1;
geom_err_dist_int = zeros(M,N);


%% Loop
% For each stalk, take the difference between the original interior data
% and the current geometric case (ellipse + some PCs).

for i = 1:m
    percent_done = 100*i/m
    
    starting_ind = (i-1)*npts + 1;    
    
    real_ext = ALL_R_ext(i,:);
    ellipse_ext = ALL_ELLIPSE_R_ext(i,:);
    real_int = ALL_R_int(i,:);
    ellipse_int = normintV2(real_ext,ALL_ELLIPSE_T(i,:),ALL_AVG_RIND_T(i));
%     ellipse_int = ALL_ELLIPSE_R_int(i,:);
    
    % Calculate the differences for the ellipse-only case
    geom_err_ext = zeros(size(real_ext'));
    geom_err_int = zeros(size(real_int'));
    for j = 1:length(geom_err_ext)
        geom_err_ext(j) = 100*(ellipse_ext(j) - real_ext(j))/(ALL_B(i)/2);
        geom_err_int(j) = 100*(ellipse_int(j) - real_int(j))/(ALL_B(i)/2);
    end
    
    % Place the ellipse errors in the array
    geom_err_dist_ext(starting_ind:(starting_ind+npts-1),1) = geom_err_ext;
    geom_err_dist_int(starting_ind:(starting_ind+npts-1),1) = geom_err_int;
    
    for j = 1:N-1
        % Calculate the geometry of the current PC case
        PC_ext = zeros(1,size(ext_rhoPCAs,1));
        for k = 1:j
            % Add all PCs up to the current PC to the ellipse in polar coordinates
            PC_ext = PC_ext + ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
        end
        
        R_ext = ALL_ELLIPSE_R_ext(i,:) - PC_ext;
        R_int = normintV2(R_ext,ALL_ELLIPSE_T(i,:),ALL_AVG_RIND_T(i));
        
        % Calculate the differences for the current PC case
        geom_err_ext = zeros(size(real_ext'));
        geom_err_int = zeros(size(real_int'));
        for k = 1:length(geom_err_ext)
            geom_err_ext(k) = 100*(R_ext(k) - real_ext(k))/ALL_B(i); % Should this be signed error or magnitude only?
            geom_err_int(k) = 100*(R_int(k) - real_int(k))/ALL_B(i);
        end
        
        % Insert the calculated differences in the error array
        geom_err_dist_ext(starting_ind:(starting_ind+npts-1),j+1) = geom_err_ext;
        geom_err_dist_int(starting_ind:(starting_ind+npts-1),j+1) = geom_err_int;

    end
    
end

    

end