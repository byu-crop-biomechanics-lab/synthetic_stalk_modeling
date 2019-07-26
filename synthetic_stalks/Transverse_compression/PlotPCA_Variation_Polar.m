function [] = PlotPCA_Variation_Polar(Theta, RprinComps_ext, Rcoeffs_ext, R_ext, Rexplained_ext, Nviz, scaling)

%--------------------------------------------------------------
% FILE: PlotPCA_Variation_Polar.m
% AUTHOR: Ryan Larson
% DATE: 5/3/19
% 
% Based on Dr. Cook's PlotPCA_Variation.m
% 
% PURPOSE: Visualization of principal components in polar coordinates
% 
% INPUTS: 
%           RprinComps - the principle components matrix from a polar data set
%           Theta - an array of theta values with the same number of rows
%           as there are PCs (the values should be the same in all rows)
%           Rcoeffs - the coefficients from PCA (referred to as "scores" in the MATLAB documentation)
%           Rmu - a vector of means from PCA
%           Nviz - an index vector, indicating the PCs to be plotted
%           scaling - an plotting option:
%                       0 - (default) uses the median coefficient value to scale each PC. This shows the relative influence of each PC 
%                       'equal', scales each PC equally. This is useful to visualize the shape itself, while the default shows the influence of each PC
%
% VERSION HISTORY
% V1 - 
% V2 - 
% V3 - 
% 
%--------------------------------------------------------------
close;

%% Exterior NEPCs

N = size(RprinComps_ext,1);             % the size of rho vectors
Means = R_ext;
Amps = median(abs(Rcoeffs_ext),1);      % Amplitudes of variation
A = 1;

% RprinComps = Rcoeffs*RprinComps';
% RprinComps = RprinComps';
avgRcoeffs = zeros(1,size(Rcoeffs_ext,2));
for i = 1:length(avgRcoeffs)
    avgRcoeffs(i) = mean(Rcoeffs_ext(:,i));
end

figure('units','normalized','outerposition',[0 0 1 1])

for i = Nviz
    
    if strcmp(scaling,'equal')
        Amps = Amps*0 + 1;
        A = 5;
        Amp = max(avgRcoeffs(i)*RprinComps_ext(:,i));
        RprinComps_ext(:,i) = avgRcoeffs(i)*RprinComps_ext(:,i)/Amp;
    end
    
    for t = linspace(0,4*pi,400)
        polarscatter(Theta, Means,'.','MarkerFaceColor',[1 1 1]*0.0)
        hold on
        polarscatter(Theta, Means - A*Amps(i)*sin(t)*RprinComps_ext(:,i)','o',...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', [0 0 0], ...
            'MarkerFaceAlpha', 0.3);
        hold on
        ax = gca;
        ax.RLim = [0 15];
        
        NEPC_No = sprintf('NEPC No. %d',i);
        percent_explained = sprintf('%0.2f%% Non-ellipse variation',Rexplained_ext(i));
        title([{NEPC_No,percent_explained}]);
        hold off
        pause(0.01)
        
    end
    
end
