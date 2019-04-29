load five_sections_R.mat

[extPCAs, extcoeffs, extPCA_variances, exttstat, extexplained, extvarMeans] = pca(DATA_ext);
[intPCAs, intcoeffs, intPCA_variances, inttstat, intexplained, intvarMeans] = pca(DATA_int);

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
