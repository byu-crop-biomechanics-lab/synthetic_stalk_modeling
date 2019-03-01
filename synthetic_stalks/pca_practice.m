%% pca_practice.m
% RL - 2/22/2019
% Script for calculating the principal components of a data set generated
% by stalk_cross_sections.m, above a chosen threshold.
clear;

% Make sure to run stalk_cross_sections.m before running this script
load cross_sections.mat

% Choose a threshold value below which to reject principal components that
% don't contribute enough to the overall shape of the data.
threshold = 4;

% Organize data from cross_sections.mat
X = sections(:,:,1);
Y = sections(:,:,2);

% Run PCA analysis on the x and y data
[coeffx, scorex, latentx, ~, explainedx, mux] = pca(X,'Centered',true,'VariableWeights','variance');
[coeffy, scorey, latenty, ~, explainedy, muy] = pca(Y,'Centered',true,'VariableWeights','variance');
% [coeffx, scorex, latentx, ~, explainedx, mux] = pca(X);
% [coeffy, scorey, latenty, ~, explainedy, muy] = pca(Y);

% Bar plots of normalized latent vectors for visualization of relative
% contribution of principal components (there are N principal components,
% but only the first few have much influence, so only the first 10 are
% plotted)
% figure(1);
% bar(explainedx);
% title('Principal components in x');
% xlim([0,10]);
% xlabel('Principal Component');
% ylabel('% Variance explained by PC');
% 
% figure(2);
% bar(explainedy);
% title('Principal components in y');
% xlim([0,10]);
% xlabel('Principal Component');
% ylabel('% Variance explained by PC');

% Count up the number of principal components that account for more than
% the percentage called out by threshold:
countx = 0;
county = 0;
for i = 1:length(explainedx)
    if explainedx(i) > threshold
        countx = countx + 1;
    end
    if explainedy(i) > threshold
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
% plot(sumx,'LineWidth',2);
title('Principal components from x data');
hold off

subplot(1,3,2);
hold on
for i = 1:county
    plot(coeffy(:,i));
    sumy = sumy + coeffy(:,i);
end
% plot(sumy,'LineWidth',2);
title('Principal components from y data');
hold off

% Reconstruct the original data using only the principal components that
% matter, based on threshold value
xapprox = scorex(:,1:countx)*coeffx(:,1:countx)';
yapprox = scorey(:,1:county)*coeffy(:,1:county)';



% This plots the un-scaled cross section shape that the principal
% components generate
subplot(1,3,3);
% plot(sumx,sumy);
sec = 5;    % choose which approximated section shape to show
% plot(xapprox(sec,:),yapprox(sec,:));
plot(sumx,sumy);