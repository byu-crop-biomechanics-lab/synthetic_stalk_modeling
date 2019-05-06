function [] = PlotPCA_Variation_Polar(RprinComps, Theta, Rcoeffs, Rmu, Nviz, scaling)

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

N = size(RprinComps,1);             % the size of rho vectors
Means = Rmu;
Amps = median(abs(Rcoeffs),1);      % Amplitudes of variation
A = 1;

% RprinComps = Rcoeffs*RprinComps';
% RprinComps = RprinComps';
avgRcoeffs = zeros(1,size(Rcoeffs,1));
for i = 1:length(avgRcoeffs)
    avgRcoeffs(i) = mean(Rcoeffs(:,i));
end

figure('units','normalized','outerposition',[0 0 1 1])

for i = Nviz
    
    if strcmp(scaling,'equal')
        Amps = Amps*0 + 1;
        A = 15;
        Amp = max(avgRcoeffs(i)*RprinComps(:,i));
        RprinComps(:,i) = avgRcoeffs(i)*RprinComps(:,i)/Amp;
    end
    
    for t = linspace(0,4*pi,400)
        polarscatter(Theta(i,:)', Means,'.','MarkerFaceColor',[1 1 1]*0.0)
        hold on
        polarscatter(Theta(i,:)', Means - A*Amps(i)*sin(t)*RprinComps(:,i)','o',...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', [0 0 0], ...
            'MarkerFaceAlpha', 0.3);
        hold on
        ax = gca;
        ax.RLim = [0 200];
        
        title(['PC No. ',num2str(i)])
        hold off
        pause(0.01)
        
    end
    
end
