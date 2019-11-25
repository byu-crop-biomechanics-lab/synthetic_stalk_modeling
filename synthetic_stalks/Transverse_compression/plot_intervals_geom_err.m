function [intervals] = plot_intervals_geom_err(geom_err_dist)
% Plot the mean and two standard deviation lines for the geom_err_dist data

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

% % Basic plot
% figure(1);
% hold on
% plot(numNEPCs,intervals(1,:));
% plot(numNEPCs,intervals(2,:));
% plot(numNEPCs,intervals(3,:));
% plot(numNEPCs,intervals(4,:));
% plot(numNEPCs,intervals(5,:));
% hold off


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

figure(2);
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