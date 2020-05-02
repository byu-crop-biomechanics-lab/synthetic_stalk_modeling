function [intervals] = plot_intervals_geom_err(geom_err_dist)
% Plot the mean and two standard deviation lines for the geom_err_dist data

hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

% Create the output data structure
intervals = zeros(5,size(geom_err_dist,2));

% Create numPCs vector for plotting
numPCs = linspace(0,size(geom_err_dist,2)-1,size(geom_err_dist,2));


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
% plot(numPCs,intervals(1,:));
% plot(numPCs,intervals(2,:));
% plot(numPCs,intervals(3,:));
% plot(numPCs,intervals(4,:));
% plot(numPCs,intervals(5,:));
% hold off


%% Patch plots for paper
% Get x values (same for all patches)
flippedPCs = fliplr(numPCs);
patchPCs = [numPCs flippedPCs];
for i = 1:length(patchPCs)
    patchPCs(i) = log10(patchPCs(i))/log10(5);
end

for i = 1:length(numPCs)
    numPCs(i) = log10(numPCs(i))/log10(5);
end
% patchPCs = log(5,patchPCs);

ellipse_plotval = -0.6;

patchPCs(1) = ellipse_plotval;
patchPCs(end) = ellipse_plotval;
numPCs(1) = ellipse_plotval;

%% Log5 plot (excludes data for ellipse approximation)
% 2 standard deviations patch
flipped_two_stdevs = fliplr(intervals(5,:));
two_stdevs = [intervals(1,:) flipped_two_stdevs];

% 1 standard deviation patch
flipped_one_stdev = fliplr(intervals(4,:));
one_stdev = [intervals(2,:) flipped_one_stdev];

% Ellipse case scatter data
% ellipse_scatter = intervals(:,1);
% ellipse_x_scatter = zeros(size(ellipse_scatter));

figure(1);
hold on
patch(patchPCs,two_stdevs,[0.75 0.75 0.75],'EdgeColor','none');
patch(patchPCs,one_stdev,[0.5 0.5 0.5],'EdgeColor','none');
plot(numPCs,intervals(3,:),'LineWidth',2,'Color','k');
% scatter(ellipse_x_scatter,ellipse_scatter,'k','filled');
hold off
xlabel('Log5(PCs Included)');
xlim([ellipse_plotval log10(360)/log10(5)]);
ylabel('Geometric Error'); % Change ylabel format to include % signs
ytickformat('percentage');
title('Progression of Geometric Accuracy');
legend('2\sigma','\sigma','Mean');
% set(gca,'XScale','log','YScale','linear')

% %% Linear plot (includes data from ellipse approximation)
% % Create numPCs vector for plotting
% numPCs = linspace(1,size(geom_err_dist,2),size(geom_err_dist,2));
% numPCs = numPCs - 1;
% 
% % Get x values (same for all patches)
% flippedPCs = fliplr(numPCs);
% patchPCs = [numPCs flippedPCs];
% % for i = 1:length(patchPCs)
% %     patchPCs(i) = log10(patchPCs(i))/log10(5);
% % end
% 
% % for i = 1:length(numPCs)
% %     numPCs(i) = log10(numPCs(i))/log10(5);
% % end
% 
% % 2 standard deviations patch
% flipped_two_stdevs = fliplr(intervals(5,:));
% two_stdevs = [intervals(1,:) flipped_two_stdevs];
% 
% % 1 standard deviation patch
% flipped_one_stdev = fliplr(intervals(4,:));
% one_stdev = [intervals(2,:) flipped_one_stdev];
% 
% figure(2);
% hold on
% patch(patchPCs,two_stdevs,[0.75 0.75 0.75],'EdgeColor','none');
% patch(patchPCs,one_stdev,[0.5 0.5 0.5],'EdgeColor','none');
% plot(numPCs,intervals(3,:),'LineWidth',2,'Color','k');
% hold off
% xlabel('PCs Included');
% ylabel('Geometric Error'); % Change ylabel format to include % signs
% ytickformat('percentage');
% title('Progression of Geometric Accuracy');
% legend('2\sigma','\sigma','Mean');


end