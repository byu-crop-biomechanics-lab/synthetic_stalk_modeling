function [sensitivities,avg_sensitivities,stderror,normalized_sensitivities] = get_sensitivity_results(Data,pct_change)
% get_reaction_results.m: Take the text file structure of the NEPC
% finite-element results and reorganize the data into an array, with the
% rows being the cross-section number and the columns being the case
% number. Cases 0-6 are used.
% clear; close;
close all
load(Data,'AllResults');

% % CONVERT FROM MICROMETER SCALE TO MILLIMETER SCALE FOR LOOKING AT ACTUAL
% % VALUES
% 
% rows = 50;
% cols = 11;
% 
% % ylimits = 0;
% % lowlim = -1.6;
% % uplim = 0.1;
% 
% Results_new = NaN(rows,cols);
% 
% for i = 1:size(Results,1)
%     if i == 1
%         row = 1;
%     elseif Results(i,1) ~= Results(i-1,1)
%         row = row + 1;
%     end
% %     row = round(Results(i,1),0);
%     col = round(Results(i,2),0) + 1;
%     
%     Results_new(row,col) = abs(Results(i,4));
% end

percents = zeros(size(AllResults));

for i = 1:size(AllResults,1)
    for j = 1:size(AllResults,2)
        percents(i,j) = (AllResults(i,j)/AllResults(i,1))*100;
    end
end

sensitivities = NaN(size(AllResults));

for i = 1:size(sensitivities,1)
    for j = 1:size(sensitivities,2)
        sensitivities(i,j) = ((AllResults(i,j) - AllResults(i,1))/AllResults(i,1))*100;
    end
end

avg_sensitivities = nanmean(sensitivities);
stdev = nanstd(sensitivities);
stderror = zeros(size(stdev));
n = sum(~isnan(sensitivities),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev)
    stderror(i) = stdev(i)/sqrt(n(i));
end

stderror = stderror';
avg_sensitivities = abs(avg_sensitivities);
avg_sensitivities = avg_sensitivities';
normalized_sensitivities = avg_sensitivities/pct_change;
normalized_stderror = stderror/pct_change;

% Avg_Sensitivities = categorical(avg_sensitivities);

% Local sensitivity bar chart
figure(1);
caselabels = {'Base'; 'A'; 'B'; 'T'; 'E_r'; 'E_p'; 'NEPC 1';...
    'NEPC 2'; 'NEPC 3'; 'NEPC 4'; 'NEPC 5'};
CaseLabels = categorical(caselabels);

T = table(avg_sensitivities,CaseLabels,stderror,normalized_sensitivities,normalized_stderror);

% Sort table values by sensitivities
T = sortrows(T,'avg_sensitivities','descend');

bar(T.avg_sensitivities);
ax = gca;
ax.XTickLabel = T.CaseLabels;

% bar(T.Avg_Sensitivities,'FaceColor',[0.75,0.75,0.75]);
% set(gca,'xticklabel',T.CaseLabels);
% if ylimits == 1
%     ylim([lowlim,uplim]);
% end
title('Average Parameter Sensitivities (5% Change)');
xlabel('Case');
ylabel('Local Sensitivity (% of Base Response)');

hold on

er = errorbar(1:11,T.avg_sensitivities,T.stderror,T.stderror);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

hold off


% Normalized sensitivity bar chart
figure(2);
bar(T.normalized_sensitivities);
ax = gca;
ax.XTickLabel = T.CaseLabels;

% caselabels = {'Base'; 'A'; 'B'; 'T'; 'E_r'; 'E_p'; 'NEPC 1';...
%     'NEPC 2'; 'NEPC 3'; 'NEPC 4'; 'NEPC 5';};
% bar(normalized_sensitivities,'FaceColor',[0.75,0.75,0.75]);
% set(gca,'xticklabel',caselabels);
% if ylimits == 1
%     ylim([lowlim,uplim]);
% end
title('Normalized Average Parameter Sensitivities (5% Change)');
xlabel('Case');
ylabel('Normalized Sensitivity');

hold on

er = errorbar(1:11,T.normalized_sensitivities,T.normalized_stderror,T.normalized_stderror);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

hold off

end