%% pca_practice.m
% RL - 2/22/2019
% Script for calculating the principal components of a data set generated
% by stalk_cross_sections.m, above a chosen threshold.
% clear;
close all;

% Make sure to run stalk_cross_sections.m before running this script
load cross_sections.mat

% Choose a threshold value below which to reject principal components that
% don't contribute enough to the overall shape of the data.
threshold = 95;

% Organize data from cross_sections.mat
X = sections(:,:,1);
Y = sections(:,:,2);

% Run PCA analysis on the x and y data
[xPCAs, xcoeffs, xPCA_variances, xtstat, xexplained, xvarMeans] = pca(X);
[yPCAs, ycoeffs, yPCA_variances, ytstat, yexplained, yvarMeans] = pca(Y);

% Count up the number of principal components that account for more than
% the percentage called out by threshold:
xcount = 0;
ycount = 0;
xsumexp = 0;
ysumexp = 0;
xlast = 0;
ylast = 0;
for i = 1:length(xexplained)
    % Step through the explained vectors and add up percentages until they
    % exceed the desired threshold
    xsumexp = xsumexp + xexplained(i);
    if xsumexp < threshold 
        xcount = xcount + 1;
    elseif xsumexp > threshold && xlast < threshold
        xcount = xcount + 1;
    end
    xlast = xsumexp;
    
    ysumexp = ysumexp + yexplained(i);
    if ysumexp < threshold
        ycount = ycount + 1;
    elseif ysumexp > threshold && ylast < threshold
        ycount = ycount + 1;
    end
    ylast = ysumexp;
end

% Determine the total percentage of the data that is captured by the chosen
% principal components (will be above chosen threshold)
PCcapturex = 0;
PCcapturey = 0;
for i = 1:xcount
    PCcapturex = PCcapturex + xexplained(i);
end
for i = 1:ycount
    PCcapturey = PCcapturey + yexplained(i);
end

% Holding vectors for sums of important principal components
sumx = zeros(size(xPCAs(:,1)));
sumy = zeros(size(yPCAs(:,1)));

minangle = lowind*thetastep;
maxangle = upind*thetastep;
angles = linspace(minangle,maxangle,size(xPCAs,2));

% Plot the principal components in x and y that are above the chosen
% threshold. Thick line is the sum of the components. 
figure('Position',[75, 250, 1800, 500]);
subplot(1,2,1);
hold on
for i = 1:xcount
%     plot(xPCAs(lowind:upind,i));    % This line only works if notchPCA.m has been run before
    plot(xPCAs(:,i));
    sumx = sumx + xPCAs(:,i);
end
% plot(sumx,'LineWidth',2);
str = sprintf('PCs from x data (%0.2f%% captured)',PCcapturex);
title(str);
xlabel('Index (centered around \pi rad)');
ylabel('PC value');
legendstrx = cell(1,xcount);
for i = 1:xcount
    str = sprintf('PC%d',i);
    legendstrx{i} = str;
end
legend(legendstrx);
hold off

subplot(1,2,2);
hold on
for i = 1:ycount
%     plot(yPCAs(lowind:upind,i));    % This line only works if notchPCA.m has been run before
    plot(yPCAs(:,i));
    sumy = sumy + yPCAs(:,i);
end
% plot(sumy,'LineWidth',2);
str = sprintf('PCs from y data (%0.2f%% captured)',PCcapturey);
title(str);
xlabel('Index (centered around \pi rad)');
ylabel('PC value');
legendstry = cell(1,ycount);
for i = 1:ycount
    str = sprintf('PC%d',i);
    legendstry{i} = str;
end
legend(legendstry);
hold off

% Reconstruct the original data using only the principal components that
% matter, based on threshold value
xapprox = xcoeffs(:,1:xcount)*xPCAs(:,1:xcount)';
yapprox = ycoeffs(:,1:ycount)*yPCAs(:,1:ycount)';

% This plots the un-scaled cross section shape that the principal
% components generate
% subplot(1,,3);
% plot(sumx(lowind:upind),sumy(lowind:upind));
% % plot(sumx,sumy);
% title('Sum of principal components in x and y');

% sec = 10;    % choose which approximated section shape to show
% plot(xapprox(sec,:),yapprox(sec,:));
% title('Reconstructed data from PCs');
