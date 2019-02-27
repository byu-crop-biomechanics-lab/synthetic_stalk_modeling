%% pca_practice.m
% RL - 2/22/2019
% Script for calculating the principal components of a data set generated
% by stalk_cross_sections.m, above a chosen threshold.
clear;

% Make sure to run stalk_cross_sections.m before running this script
load cross_sections.mat

% Choose a threshold value below which to reject principal components that
% don't contribute enough to the overall shape of the data.
threshold = 0.1;

% Organize data from cross_sections.mat
X = sections(:,:,1);
Y = sections(:,:,2);

% Run PCA analysis on the x and y data
[coeffx, scorex, latentx] = pca(X);
[coeffy, scorey, latenty] = pca(Y);

% Create vectors that normalize latentx and latenty
pctlatx = zeros(size(latentx));
pctlaty = zeros(size(latenty));

for i = 1:length(pctlatx)
    pctlatx(i) = latentx(i)/sum(latentx);
end

for i = 1:length(pctlaty)
    pctlaty(i) = latenty(i)/sum(latenty);
end

% Bar plots of normalized latent vectors for visualization of relative
% contribution of principal components (there are N principal components,
% but only the first few have much influence, so only the first 10 are
% plotted)
figure(1);
bar(pctlatx);
title('Principal components in x');
xlim([0,10]);

figure(2);
bar(pctlaty);
title('Principal components in y');
xlim([0,10]);

% if length(pctlatx) ~= length(pctlaty)
%     error('latent vectors are not the same length');
% end

% Count up the number of principal components that account for more than
% the percentage called out by threshold:
countx = 0;
county = 0;
for i = 1:length(pctlatx)
    if pctlatx(i) > threshold
        countx = countx + 1;
    end
    if pctlaty(i) > threshold
        county = county + 1;
    end
end

% Holding vectors for sums of important principal components
sumx = zeros(size(coeffx(:,1)));
sumy = zeros(size(coeffy(:,1)));

% Plot the principal components in x and y that are above the chosen
% threshold. Thick line is the sum of the components. 
figure(3);
subplot(1,3,1);
hold on
for i = 1:countx
    plot(coeffx(:,i));
    sumx = sumx + coeffx(:,i);
end
plot(sumx,'LineWidth',2);
title('Principal components from x data');
hold off

subplot(1,3,2);
hold on
for i = 1:county
    plot(coeffy(:,i));
    sumy = sumy + coeffy(:,i);
end
plot(sumy,'LineWidth',2);
title('Principal components from y data');
hold off

% This plots the un-scaled cross section shape that the principal
% components generate
subplot(1,3,3);
plot(sumx,sumy);