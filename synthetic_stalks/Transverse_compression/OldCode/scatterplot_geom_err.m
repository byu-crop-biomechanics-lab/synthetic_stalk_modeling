function scatterplot_geom_err(geom_err_dist)
% Create a scatter plot of the output from get_geom_err_dist.m.

close all

numNEPCs = zeros(size(geom_err_dist));
for i = 1:size(numNEPCs,1)
    numNEPCs(i,:) = linspace(0,size(numNEPCs,2)-1,size(numNEPCs,2));
end

% Scatter plot of total errors
hold on
for i = 1:size(numNEPCs,1)
    scatter(numNEPCs(i,:),geom_err_dist(i,:));    
end
hold off

end