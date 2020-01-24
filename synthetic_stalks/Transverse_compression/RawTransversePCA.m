function RawTransversePCA(AllSlicesPCA)
% FILENAME: RawTransversePCA.m
% AUTHOR: Ryan Larson
% DATE: 1/24/2020
%
% PURPOSE: 
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
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load(AllSlicesPCA);

[ext_rhoPCAs, ext_rhocoeffs, ext_rhoPCA_variances, ext_rhotstat, ext_rhoexplained, ext_rhovarMeans] = pca(ALL_R_ext,'Centered',false);

polarplot(ELLIPSE_T(1,:),ext_rhoPCAs(:,1));
hold on
polarplot(ELLIPSE_T(1,:),ext_rhoPCAs(:,2));
polarplot(ELLIPSE_T(1,:),ext_rhoPCAs(:,3));
polarplot(ELLIPSE_T(1,:),ext_rhoPCAs(:,4));
polarplot(ELLIPSE_T(1,:),ext_rhoPCAs(:,5));

hold off
legend('PC1','PC2','PC3','PC4','PC5');




end