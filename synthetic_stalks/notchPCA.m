close all

load cross_sections.mat

X = sections(:,:,1);
Y = sections(:,:,2);
Xnotch = X;
Ynotch = Y;

N = size(X,2);
threshold = 90;

% Flip data so notch is always on the right side
X = -X;
Y = -Y;

thetastep = (2*pi)/N;

notchwidths = zeros(size(X,1),1);
notchlocs = zeros(size(X,1),2);
% notchranges = notchwidths;
% notchrangelocs = notchlocs;

for i = 1:length(notchwidths)
    sel = (max(X(i,:))-min(X(i,:)))/200;
%     peakfinder(X(i,:),sel);
%     pause();
%     close;
    [peakloc,~] = peakfinder(X(i,:),sel);
%     if length(peakloc) ~= 2
%         length(peakloc)
%         error('peakfinder function didn''t find two peaks');
%     end
    if length(peakloc) ~= 2
        notchlocs(i,1) = NaN;
        notchlocs(i,2) = NaN;
        notchwidths(i) = NaN;
    else
    notchlocs(i,1) = peakloc(1);
    notchlocs(i,2) = peakloc(2);
    notchrange = notchlocs(i,2) - notchlocs(i,1);
%     % Adjust indices of notch locations by half the range
%     notchlocs(i,1) = notchlocs(i,1) - ceil(notchrange/2);
%     notchlocs(i,2) = notchlocs(i,2) + ceil(notchrange/2);
    notchwidths(i) = (notchlocs(i,2) - notchlocs(i,1))*thetastep;
    end
end

% histogram(notchwidths,15)

toploc = max(notchlocs(:,2));
botloc = min(notchlocs(:,1));

% Cut out data on each cross section that doesn't fall between toploc and
% botloc indices
for i = 1:size(X,1)
    for j = 1:N
        if j < botloc || j > toploc
            Xnotch(i,j) = nan;
            Ynotch(i,j) = nan;
        end
    end
end

% Get rid of nan values in Xnotch and Ynotch
Xnotch = rmmissing(Xnotch,2);
Ynotch = rmmissing(Ynotch,2);

% for i = 1:size(Xnotch,1)
%     plot(Xnotch(i,:),Ynotch(i,:));
%     pause();
% end

% Run PCA analysis on the x and y data
[xPCAs, xcoeffs, xPCA_variances, xtstat, xexplained, xvarMeans] = pca(Xnotch);
[yPCAs, ycoeffs, yPCA_variances, ytstat, yexplained, yvarMeans] = pca(Ynotch);

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

% Plot the principal components in x and y that are above the chosen
% threshold. Thick line is the sum of the components. 
figure('Position',[75, 250, 1800, 500]);
subplot(1,3,1);
hold on
for i = 1:xcount
    plot(xPCAs(:,i));
    sumx = sumx + xPCAs(:,i);
end
% plot(sumx,'LineWidth',2);
str = sprintf('PCs from x data (%0.2f%% captured)',PCcapturex);
title(str);
hold off

subplot(1,3,2);
hold on
for i = 1:ycount
    plot(yPCAs(:,i));
    sumy = sumy + yPCAs(:,i);
end
% plot(sumy,'LineWidth',2);
str = sprintf('PCs from y data (%0.2f%% captured)',PCcapturey);
title(str);
hold off

% Reconstruct the original data using only the principal components that
% matter, based on threshold value
xapprox = xcoeffs(:,1:xcount)*xPCAs(:,1:xcount)';
yapprox = ycoeffs(:,1:ycount)*yPCAs(:,1:ycount)';

% This plots the un-scaled cross section shape that the principal
% components generate
subplot(1,3,3);
plot(sumx,sumy);
title('Sum of principal components in x and y');

% sec = 10;    % choose which approximated section shape to show
% plot(xapprox(sec,:),yapprox(sec,:));
% title('Reconstructed data from PCs');
