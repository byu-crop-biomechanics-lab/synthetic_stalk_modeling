function [] = PlotPCA_Variation(prinComps, coeffs, mu, Nviz, scaling)

%--------------------------------------------------------------
% FILE: PlotPCA_Variation.m
% AUTHOR: DC
% DATE: 3/8/19
% 
% PURPOSE: Visualization of prinComps
% 
% INPUTS: 
%           prinComps - the principle components matrix from a XY data set
%           coeffs - the coefficients from PCA (referred to as "scores" in the MATLAB documentation)
%           mu - a vector of means from PCA
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


N = size(prinComps,2)/2;                    % the size of individual x and y vectors
xprinComps = prinComps(1:N,:);          % x and y prinComps
yprinComps = prinComps(N+1:2*N,:);        
xMeans = mu(1:N);                       % x and y mean values   
yMeans = mu(N+1:2*N);   
Amps = median(abs(coeffs),1); %*2+0.3;               % Amplitudes of variation
A = 1;

figure('units','normalized','outerposition',[0 0 1 1])


for i = Nviz
    
    if strcmp(scaling,'equal')
        Amps = Amps*0 + 1;
        A = 0.1;
        Amp = max(sqrt(xprinComps(:,i).^2 + xprinComps(:,i).^2));
        xprinComps(:,i) = xprinComps(:,i)/Amp;
        yprinComps(:,i) = yprinComps(:,i)/Amp;
    end
    
    
    for t = linspace(0,4*pi,400)

        scatter(xMeans, yMeans,'.','MarkerFaceColor',[1 1 1]*0.0)
        hold on
        h = scatter(xMeans + A*Amps(i)*sin(t)*xprinComps(:,i)',yMeans + A*Amps(i)*sin(t)*yprinComps(:,i)','o',...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', [0 0 0], ...
            'MarkerFaceAlpha', 0.3);
        
        title(['PC No. ',num2str(i)])
        hold off
        xlim([-1 1]*0.8)
        ylim([-1 1]*0.6)
        if strcmp(scaling,'equal'), text(-0.25,0,['Note: Plotting is Exaggerated to Show Detail.']), end
        axis equal
        pause(0.01)
        
    end
    
end
