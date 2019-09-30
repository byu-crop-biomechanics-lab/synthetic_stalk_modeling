function [intervals] = plot_intervals_geom_err(geom_err_dist)
% Plot the mean and two standard deviation lines for the geom_err_dist data

hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

% Create the output data structure
intervals = zeros(5,size(geom_err_dist,2));

% Create numNEPCs vector for plotting
numNEPCs = linspace(0,size(geom_err_dist,2)-1,size(geom_err_dist,2));


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

% Basic plot
figure(1);
hold on
plot(numNEPCs,intervals(1,:));
plot(numNEPCs,intervals(2,:));
plot(numNEPCs,intervals(3,:));
plot(numNEPCs,intervals(4,:));
plot(numNEPCs,intervals(5,:));
hold off


%% Patch plot for paper
% Get x values (same for all patches)
flippedNEPCs = fliplr(numNEPCs);
patchNEPCs = [numNEPCs flippedNEPCs];

% 2 standard deviations patch
flipped_two_stdevs = fliplr(intervals(5,:));
two_stdevs = [intervals(1,:) flipped_two_stdevs];

% 1 standard deviation patch
flipped_one_stdev = fliplr(intervals(4,:));
one_stdev = [intervals(2,:) flipped_one_stdev];

figure(2);
hold on
patch(patchNEPCs,two_stdevs,[0.75 0.75 0.75],'EdgeColor','none');
patch(patchNEPCs,one_stdev,[0.5 0.5 0.5],'EdgeColor','none');
plot(numNEPCs,intervals(3,:),'LineWidth',2,'Color','k');
hold off
xlabel('Number of NEPCs');
ylabel('Total Geometric Error (mm)');
title('Progression of Geometric Accuracy');
legend('2\sigma','\sigma','Mean');


end