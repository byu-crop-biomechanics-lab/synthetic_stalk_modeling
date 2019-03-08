function [] = PlotPCA_Fitting(DAT, prinComps, coeffs, mu, CrossSections, Ncomps)

%--------------------------------------------------------------
% FILE: PlotPCA_Fitting.m
% AUTHOR: DC
% DATE: 3/8/19
% 
% PURPOSE: Visualization of how PCs combine to re-create original data
% 
% INPUTS: 
%           DAT - the data passed into PCA, assumed to be a combined XY data set
%           prinComps - the principle components matrix from a XY data set
%           coeffs - the coefficients from PCA (referred to as "scores" in the MATLAB documentation)
%           mu - a vector of means from PCA
%           CrossSections - which cross-sections to plot
%           Ncomps - the total number of PCs to show for each cross-section
%
% VERSION HISTORY
% V1 - 
% V2 - 
% V3 - 
% 
%--------------------------------------------------------------

N = size(DAT,2)/2;          % the number of variables in x (or y)                   
xPCAs = prinComps(1:N,:);            % the x-coordinate PCs
yPCAs = prinComps(N+1:2*N,:);        % the y-coordinate PCs
xMeans = mu(1:N);                    % the x-means
yMeans = mu(N+1:2*N);                % y-means
oldRsq = 0;

for j = CrossSections
    for i = 1:Ncomps

        xCS = DAT(j,[1:N,1]);               % x cross-section data
        yCS = DAT(j,[N+1:2*N,N+1]);
        plot(xCS, yCS, 'Color',[1 1 1]*0.7)
        hold on
        
        xPred = coeffs(j,1:i)*xPCAs([1:end,1],1:i)' + xMeans([1:end,1]);     % Predicted x-data
        yPred = coeffs(j,1:i)*yPCAs([1:end,1],1:i)' + yMeans([1:end,1]);
        plot( xPred, yPred, '.')
        
        % Calculating Rsq based on plots of CS and Pred, see below (x and y residuals are handled independently
        CS = [xCS, yCS];            % Cross-section data as a vector
        Pred = [xPred, yPred];      % Predicted data as a vector
        residuals = Pred - CS;      % residuals
        SST = sum((CS - mean(CS)).^2);  % sums squared total
        SSR = sum(residuals.^2);        % sums squared of the residuals
        Rsq = 1 - SSR/SST;              % R-squared value
        
        % Calculating Rsq based on cartesian distance (This turns out to be identical to the above!)
        SqResiduals = (xPred - xCS).^2 + (yPred - yCS).^2;     
        SSTsq = sum(xCS.^2 + yCS.^2);
        SSRsq = sum(SqResiduals);
        RsqSq = 1 - SSRsq/SSTsq;
        [maxError, maxdex] = max(sqrt(SqResiduals));  % the maximum cartesian discrepancy between real and predicted data
        r = sqrt(mean(xCS.^2 + yCS.^2));
        
        
        title([num2str(i), ' PCs'])
        text(-0.2,0.05,['Rsq = ',num2str(Rsq)])
        text(-0.2,-0.05,['Improvement: ',num2str(Rsq-oldRsq)]);
        text(-0.2,-0.15,['Max Error: ', num2str(maxError/r*100), '%'])
        
        if i == Ncomps
            text(-0.2, 0.15, ['Last Prediction.'], 'FontSize', 14)
            pause(0.5)
            
        end
        
        oldRsq = Rsq;
        pause
        hold off  
        
    end
end

