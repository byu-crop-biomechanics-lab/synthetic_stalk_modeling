close all
load cross_sections.mat

% Construct empty array that is big enough to fit all data in one panel
datasize = size(sections);
data = zeros(datasize(2),datasize(1));

% for each 
for i = 1:datasize(1)
    
    % odd/even indicator
    j = mod(i,2);
    
    % transform i to x and y indices in sections array
    ind = floor(i/2) + j;
    % odd indices get filled with x data, even indices get y data
    if j == 1
        data(:,i) = transpose(sections(ind,:,1));
    else
        data(:,i) = transpose(sections(ind,:,2));
    end
end

% need to separate this out maybe?
[PCAs,coeffs,PCA_variances,tstat,explained,varMeans] = pca(data(:,1:2));

threshold = 90;

count = 0;
sumexp = 0;
last = 0;
for i = 1:length(explained)
    sumexp = sumexp + explained(i);
    if sumexp < threshold
        count = count + 1;
    elseif sumexp > threshold && last < threshold
        count = count + 1;
    end
    last = sumexp;
end

PCcapture = 0;
for i = 1:count
    PCcapture = PCcapture + explained(i);
end

figure(1);
hold on
for i = 1:count
    plot(PCAs(:,i));
end
hold off