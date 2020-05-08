function [intervals] = plot_intervals_geom_err(geom_err_dist)
% Plot the mean and two standard deviation lines for the geom_err_dist data

% FILENAME: get_geom_err_distV2.m
% AUTHOR: Ryan Larson
% DATE: 2019
%
% PURPOSE: Gather geometric error data for each cross-section under
% examination. The inputs have been changed since version 1 to simplify and
% make the outputs more comprehensive, instead of limiting them to a single
% slice location at a time.
% 
% 
% INPUTS:
%       geom_err_dist: An array of geometric errors,  output from 
%       get_geom_err_distV2.m or get_geom_err_dist.m.
%       
% OUTPUTS:
%       intervals: 
% 
% NOTES: 
%       
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
% V1 - 
% V2 - 
% V3 - 
% -------------------------------------------------------------------------

hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

% Create the output data structure
intervals = zeros(5,size(geom_err_dist,2));

% Create numNEPCs vector for plotting
numNEPCs = linspace(1,size(geom_err_dist,2)-1,size(geom_err_dist,2)-1);


% Calculate the mean and standard deviation intervals for each geometry
% case
for i = 1:size(geom_err_dist,2)
    avg = mean(geom_err_dist(:,i));
    stdev = std(geom_err_dist(:,i));
    
    intervals(1,i) = avg + 2*stdev;
    intervals(2,i) = avg + stdev;
    intervals(3,i) = avg;
    intervals(4,i) = avg - stdev;
    intervals(5,i) = avg - 2*stdev;    
    
end


%% Patch plot for paper
% Get x values (same for all patches)
flippedNEPCs = fliplr(numNEPCs);
patchNEPCs = [numNEPCs flippedNEPCs];
for i = 1:length(patchNEPCs)
    patchNEPCs(i) = log10(patchNEPCs(i))/log10(5);
end

for i = 1:length(numNEPCs)
    numNEPCs(i) = log10(numNEPCs(i))/log10(5);
end
% patchNEPCs = log(5,patchNEPCs);

% 2 standard deviations patch
flipped_two_stdevs = fliplr(intervals(5,2:end));
two_stdevs = [intervals(1,2:end) flipped_two_stdevs];

% 1 standard deviation patch
flipped_one_stdev = fliplr(intervals(4,2:end));
one_stdev = [intervals(2,2:end) flipped_one_stdev];

figure(1);
hold on
patch(patchNEPCs,two_stdevs,[0.75 0.75 0.75],'EdgeColor','none');
patch(patchNEPCs,one_stdev,[0.5 0.5 0.5],'EdgeColor','none');
plot(numNEPCs,intervals(3,2:end),'LineWidth',2,'Color','k');
hold off
xlabel('Log5(Number of NEPCs)');
ylabel('% Geometric Error (scaled by minor diameter)');
title('Progression of Geometric Accuracy');
legend('2\sigma','\sigma','Mean');
% set(gca,'XScale','log','YScale','linear')


end