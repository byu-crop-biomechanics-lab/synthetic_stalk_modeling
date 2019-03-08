%% PCAsimple.m
% Script for examining how clearly the principal components of stalk
% features can be extracted, given prior visual extraction (i.e. notch
% alone, etc.)
clear;

nvars = 100;     % number of variables (points along shape)
nobs = 500;      % number of observations

DATA = zeros(nobs,nvars);

x = linspace(-2,2,nvars);

for i = 1:nobs
    shapecoeff = unifrnd(-0.2,0.2);
    for j = 1:nvars
        noise = unifrnd(-0.05,0.05);
        DATA(i,j) = (1 + shapecoeff)*x(j)^2 + noise;
    end
end

figure(1);
hold on
for i = 1:nobs
    plot(DATA(i,:));
%     pause(0.01);
end
hold off

[PCAs,coeffs,PCA_variances,tstat,explained,mu] = pca(DATA);

figure(2);
plot(PCAs(:,1))
hold on
plot(PCAs(:,2))