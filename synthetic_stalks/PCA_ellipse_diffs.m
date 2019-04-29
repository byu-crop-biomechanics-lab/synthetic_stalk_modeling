load DATA.mat
load ellDATA.mat

% Remove the last column in DATA_ext and DATA_int (repeats of the first
% point)
DATA_ext = DATA_ext(:,1:(end-1));
DATA_int = DATA_int(:,1:(end-1));

% Subtract the real data from the ellipse fit for the interior and exterior
ext_diff = ellDATA_ext - DATA_ext;
int_diff = ellDATA_int - DATA_int;

[extPCAs, extcoeffs, extPCA_variances, exttstat, extexplained, extvarMeans] = pca(ext_diff);
[intPCAs, intcoeffs, intPCA_variances, inttstat, intexplained, intvarMeans] = pca(int_diff);

threshold = 95;

figure(1);
bar(extexplained);
title('Exterior PCA Results');
xlabel('Principal Component');
ylabel('Percentage Explained');

figure(2);
bar(intexplained);
title('Interior PCA Results');
xlabel('Principal Component');
ylabel('Percentage Explained');

figure(3);
plot(extPCAs);
title('Exterior Principal Components');
legend('PC1','PC2','PC3','PC4');

figure(4);
plot(intPCAs);
title('Interior Principal Components');
legend('PC1','PC2','PC3','PC4');
