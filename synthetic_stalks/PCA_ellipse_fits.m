function PCA_ellipse_fits(FileName,SaveName)

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'DIFF_ext_R','DIFF_int_R','ELLIPSE_T');

% Perform PCA. 'Centered' option must be false to allow for reverse
% engineering of the original data
[ext_rhoPCAs, ext_rhocoeffs, ext_rhoPCA_variances, ext_rhotstat, ext_rhoexplained, ext_rhovarMeans] = pca(DIFF_ext_R,'Centered',false);
[int_rhoPCAs, int_rhocoeffs, int_rhoPCA_variances, int_rhotstat, int_rhoexplained, int_rhovarMeans] = pca(DIFF_int_R,'Centered',false);

ext_rhoexplained_tot = zeros(5);
int_rhoexplained_tot = zeros(5);
for i = 1:length(ext_rhoexplained_tot)
    ext_rhoexplained_tot(i) = sum(ext_rhoexplained(1:i));
    int_rhoexplained_tot(i) = sum(int_rhoexplained(1:i));
end

% figure(1);
% plot(ext_rhoexplained_tot(:,1),'-*');
% title('Exterior Non-Ellipse PCs');
% xlabel('# of PCs');
% ylabel('% Variance Explained');
% 
% figure(2);
% plot(int_rhoexplained_tot(:,1),'-*');
% title('Interior Non-Ellipse PCs');
% xlabel('# of PCs');
% ylabel('% Variance Explained');
% 
% ELLIPSE_T = ELLIPSE_T';
% theta = ELLIPSE_T(:,1);
% 
% figure(3);
% polarplot(theta,ext_rhoPCAs(:,1));
% hold on
% polarplot(theta,ext_rhoPCAs(:,2));
% polarplot(theta,ext_rhoPCAs(:,3));
% polarplot(theta,ext_rhoPCAs(:,4));
% polarplot(theta,ext_rhoPCAs(:,5));
% title('Exterior Rho Principal Components');
% legend('PC1','PC2','PC3','PC4','PC5');
% 
% figure(4);
% polarplot(theta,int_rhoPCAs(:,1));
% hold on
% polarplot(theta,int_rhoPCAs(:,2));
% polarplot(theta,int_rhoPCAs(:,3));
% polarplot(theta,int_rhoPCAs(:,4));
% polarplot(theta,int_rhoPCAs(:,5));
% title('Interior Rho Principal Components');
% legend('PC1','PC2','PC3','PC4','PC5');




% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'ELLIPSE_T','ext_rhocoeffs','ext_rhoPCAs','ext_rhoexplained',...
    'ext_rhovarMeans','int_rhoPCAs','int_rhocoeffs','int_rhoexplained','int_rhovarMeans');


end