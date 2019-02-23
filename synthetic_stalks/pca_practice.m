clear;
load cross_sections.mat

X = sections(:,:,1);
Y = sections(:,:,2);

[coeffx, scorex, latentx] = pca(X);
[coeffy, scorey, latenty] = pca(Y);

pctlatx = zeros(size(latentx));
pctlaty = zeros(size(latenty));

for i = 1:length(pctlatx)
    pctlatx(i) = latentx(i)/sum(latentx);
end

for i = 1:length(pctlaty)
    pctlaty(i) = latenty(i)/sum(latenty);
end

% figure(1);
% bar(pctlatx);
% title('Principal components in x');
% figure(2);
% bar(pctlaty);
% title('Principal components in y');

countx = 0;
county = 0;

if length(pctlatx) ~= length(pctlaty)
    error('latent vectors are not the same length');
end

threshold = 0.1;

% Count up the number of principal components that account for more than
% the percentage called out by threshold:
for i = 1:length(pctlatx)
    if pctlatx(i) > threshold
        countx = countx + 1;
    end
    if pctlaty(i) > threshold
        county = county + 1;
    end
end

sumx = zeros(size(coeffx(:,1)));
sumy = zeros(size(coeffy(:,1)));

% if countx > county
%     county = countx;
% end
% if county > countx
%     countx = county;
% end

figure(1);
subplot(1,3,1);
hold on
for i = 1:countx
    plot(coeffx(:,i));
    sumx = sumx + coeffx(:,i);
end
plot(sumx,'LineWidth',2);
title('Principal components from x data');
hold off

% figure(2);
subplot(1,3,2);
hold on
for i = 1:county
    plot(coeffy(:,i));
    sumy = sumy + coeffy(:,i);
end
plot(sumy,'LineWidth',2);
title('Principal components from y data');
hold off

% figure(3);
subplot(1,3,3);
plot(sumx,sumy);

% figure(1);
% bar(latentx);
% figure(2);
% bar(latenty);