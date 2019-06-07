% get_reaction_results.m: Take the text file structure of the NEPC
% finite-element results and reorganize the data into an array, with the
% rows being the cross-section number and the columns being the case
% number. Cases 0-6 are used.
clear; close;

load NEPC_results.mat

rows = 50;
cols = 7;

Results_new = NaN(rows,cols);

for i = 1:size(Results,1)
    row = round(Results(i,1),0);
    col = round(Results(i,2),0) + 1;
    
    Results_new(row,col) = abs(Results(i,4));
end

percents = zeros(size(Results_new));

for i = 1:size(Results_new,1)
    for j = 1:size(Results_new,2)
        percents(i,j) = (Results_new(i,j)/Results_new(i,1))*100;
    end
end

avg = nanmean(percents);
stdev = nanstd(percents);
stderror = zeros(size(stdev));
n = sum(~isnan(percents),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev)
    stderror(i) = stdev(i)/sqrt(n(i));
end

caselabels = {'Real'; 'Ellipse'; 'NEPC 1'; 'NEPC 1-2'; 'NEPC 1-3'; 'NEPC 1-4'; 'NEPC 1-5'};

bar(avg,'FaceColor',[0.75,0.75,0.75]);
set(gca,'xticklabel',caselabels);
ylim([98,100.5]);
title('NEPC Response Progression (Cumulative)');
xlabel('Case');
ylabel('Reaction Force (% of Real Response)');

hold on

er = errorbar(1:7,avg,stderror,stderror);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

% Define noisy region as +/- 0.2% of the real cross-section response
yline(99.8,':k');
yline(100.2,':k');

hold off