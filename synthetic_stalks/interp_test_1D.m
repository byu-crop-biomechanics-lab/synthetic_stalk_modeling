% interp_test_1D: Testing capabilities of interpolation functions to follow
% real data similar to that in Figure 5 of Robertson et al 2015

close all

ptCloud = pcread('exterior.ply');

% Break the point cloud into vectors
X = ptCloud.Location(:,1);
Y = ptCloud.Location(:,2);
Z = ptCloud.Location(:,3);
Zslice = [];
XYslice = []; 
countvector = [];
count = 0;
%% Restructuring data to make interp3 work the way I intend:
for i = 1:length(Z)
    if i == 1
        count = 1;
    elseif Z(i) == Z(i-1)
        count = count + 1;
    elseif Z(i) ~= Z(i-1)
        countvector = [countvector;count];
        count = 0;
    end
end

stdev = std(countvector)
avg = mean(countvector)
histogram(countvector)

% From the results above, it is clear that the real method doesn't use the
% same number of points per slice, and thus the key section interpolation
% method can't apply. Moving this exploration on to synthetic stalks
% only...

% plot3(X,Y,Z);
% hold on
% scatter(Xc,Yc,Zc,'filled');






% %% Other method
% % Load in data:
% % load 'minor_diam_fig5.mat'
% 
% % Make sure points are ordered properly...
% 
% X = minor_diam_fig5(:,1);
% Y = minor_diam_fig5(:,2);
% % N = 5;
% % P = polyfit(X,Y,N)
% % % Remember how to plot the polynomial
% % x = linspace(0,max(X));
% % fitline = P(1)*x.^5 + P(2)*x.^4 + P(3)*x.^3 + P(4)*x.^2 + P(5)*x + P(6);
% 
% Xselect = X(42:end);
% Yselect = Y(42:end);
% Xq = 405:.125:422;
% Vq = interp1(Xselect,Yselect,Xq,'pchip');
% 
% plot(Xselect,Yselect);
% hold on
% % plot(x,fitline);
% plot(Xq,Vq);
% 
