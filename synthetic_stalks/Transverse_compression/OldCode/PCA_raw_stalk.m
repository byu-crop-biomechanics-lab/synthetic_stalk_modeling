function PCA_raw_stalk(FileName)

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'ext_rhoDCSR','ext_xDCSR','ext_yDCSR','int_rhoDCSR','int_xDCSR','int_yDCSR','tDCSR','avg_rind_thickness');

% Convert exterior data to PCA-friendly form
DATA_ext_X = squeeze(ext_xDCSR);
DATA_ext_X = DATA_ext_X';
DATA_ext_Y = squeeze(ext_yDCSR);
DATA_ext_Y = DATA_ext_Y';
DATA_ext_rho = squeeze(ext_rhoDCSR);
DATA_ext_rho = DATA_ext_rho';

% Convert interior data to PCA-friendly form
DATA_int_X = squeeze(int_xDCSR);
DATA_int_X = DATA_int_X';
DATA_int_Y = squeeze(int_yDCSR);
DATA_int_Y = DATA_int_Y';
DATA_int_rho = squeeze(int_rhoDCSR);
DATA_int_rho = DATA_int_rho';

[ext_XPCAs, ext_Xcoeffs, ext_XPCA_variances, ext_Xtstat, ext_Xexplained, ext_XvarMeans] = pca(DATA_ext_X,'Centered',false);
[ext_YPCAs, ext_Ycoeffs, ext_YPCA_variances, ext_Ytstat, ext_Yexplained, ext_YvarMeans] = pca(DATA_ext_Y,'Centered',false);
[int_XPCAs, int_Xcoeffs, int_XPCA_variances, int_Xtstat, int_Xexplained, int_XvarMeans] = pca(DATA_int_X,'Centered',false);
[int_YPCAs, int_Ycoeffs, int_YPCA_variances, int_Ytstat, int_Yexplained, int_YvarMeans] = pca(DATA_int_Y,'Centered',false);
[ext_rhoPCAs, ext_rhocoeffs, ext_rhoPCA_variances, ext_rhotstat, ext_rhoexplained, ext_rhovarMeans] = pca(DATA_ext_rho,'Centered',false);
[int_rhoPCAs, int_rhocoeffs, int_rhoPCA_variances, int_rhotstat, int_rhoexplained, int_rhovarMeans] = pca(DATA_int_rho,'Centered',false);

% figure(1);
% bar(ext_Xexplained);
% title('Exterior X PCA Results');
% xlabel('Principal Component');
% ylabel('Percentage Explained');
% 
% figure(2);
% bar(ext_Yexplained);
% title('Exterior Y PCA Results');
% xlabel('Principal Component');
% ylabel('Percentage Explained');
% 
% figure(3);
% bar(int_Xexplained);
% title('Interior X PCA Results');
% xlabel('Principal Component');
% ylabel('Percentage Explained');
% 
% figure(4);
% bar(int_Yexplained);
% title('Interior Y PCA Results');
% xlabel('Principal Component');
% ylabel('Percentage Explained');

figure(5);
bar(ext_rhoexplained);
title('Exterior Rho PCA Results');
xlabel('Principal Component');
ylabel('Percentage Explained');

figure(6);
bar(int_rhoexplained);
title('Interior Rho PCA Results');
xlabel('Principal Component');
ylabel('Percentage Explained');

% figure(7);
% plot(ext_XPCAs(:,1:3));
% title('Exterior X Principal Components');
% legend('PC1','PC2','PC3');
% 
% figure(8);
% plot(ext_YPCAs(:,1:3));
% title('Exterior Y Principal Components');
% legend('PC1','PC2','PC3');
% 
% figure(9);
% plot(int_XPCAs(:,1:3));
% title('Interior X Principal Components');
% legend('PC1','PC2','PC3');
% 
% figure(10);
% plot(int_YPCAs(:,1:3));
% title('Interior Y Principal Components');
% legend('PC1','PC2','PC3');

figure(11);
polarplot(tDCSR,ext_rhoPCAs(:,1));
title('Exterior Rho Principal Components');
legend('PC1');

figure(12);
polarplot(tDCSR,int_rhoPCAs(:,1));
title('Interior Rho Principal Components');
legend('PC1');
end