function [SensitivityTable] = get_sensitivity_resultsV2(Data,pct_change)
% FILENAME: get_sensitivity_resultsV2.m
% AUTHOR: Ryan Larson
% DATE: 2/3/2020
%
% PURPOSE: Calculate sensitivities and plot them for each variable
% 
% 
% INPUTS:
%       
%       
% OUTPUTS:
%       - 
%
%
% NOTES: 
%       - 
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

close;
load(Data,'ResultsCell');

% CONVERT FROM MICROMETER SCALE TO MILLIMETER SCALE FOR LOOKING AT ACTUAL
% VALUES

rows = 70;
cols = 11;

Results_new = NaN(rows,cols);

for j = 1:length(ResultsCell)
    Results_temp = NaN(rows,cols);
    Results = cell2mat(ResultsCell(j));

    for i = 1:size(Results,1)
        if i == 1
            row = 1;
        elseif Results(i,1) ~= Results(i-1,1)
            row = row + 1;
        end
    %     row = round(Results(i,1),0);
        col = round(Results(i,2),0) + 1;

        Results_temp(row,col) = abs(Results(i,4));
    end

    if j == 1
        Results_new = Results_temp;
    else
        Results_new = [Results_new; Results_temp];
    end
end

percents = zeros(size(Results_new));

for i = 1:size(Results_new,1)
    for j = 1:size(Results_new,2)
        percents(i,j) = (Results_new(i,j)/Results_new(i,1))*100;
    end
end

sensitivities = NaN(size(Results_new));

for i = 1:size(sensitivities,1)
    for j = 1:size(sensitivities,2)
        sensitivities(i,j) = ((Results_new(i,j) - Results_new(i,1))/Results_new(i,1));
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
abs_avg_sensitivities = abs(avg_sensitivities);
abs_avg_sensitivities = abs_avg_sensitivities';
avg_sensitivities = avg_sensitivities';
normalized_sensitivities = avg_sensitivities/pct_change;
normalized_stderror = stderror/pct_change;

% Avg_Sensitivities = categorical(avg_sensitivities);

% Local sensitivity bar chart

caselabels = {'Base'; 'Major Diameter'; 'Minor Diameter'; 'Rind Thickness';...
    'Rind Modulus'; 'Pith Modulus'; 'PC 1';'PC 2'; 'PC 3'; 'PC 4'; 'PC 5'};
CaseLabels = categorical(caselabels);

SensitivityTable = table(avg_sensitivities,CaseLabels,stderror,normalized_sensitivities,normalized_stderror,abs_avg_sensitivities);

% Sort table values by sensitivities
SensitivityTable = sortrows(SensitivityTable,'abs_avg_sensitivities','descend');
SensitivityTable(end,:) = [];

% figure(1);
% bar(SensitivityTable.avg_sensitivities);
% ax = gca;
% ax.XTickLabel = SensitivityTable.CaseLabels;
% set(ax,'XTickLabelRotation',90);
% 
% % bar(T.Avg_Sensitivities,'FaceColor',[0.75,0.75,0.75]);
% % set(gca,'xticklabel',T.CaseLabels);
% % if ylimits == 1
% %     ylim([lowlim,uplim]);
% % end
% title('Average Parameter Sensitivities (10% Change)');
% xlabel('Case');
% ylabel('Local Sensitivity (% of Base Response)');
% 
% hold on
% 
% er = errorbar(1:10,SensitivityTable.avg_sensitivities,SensitivityTable.stderror,SensitivityTable.stderror);
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% er.LineWidth = 0.5;
% 
% hold off


% Normalized sensitivity bar chart
figure(2);
bar(abs(SensitivityTable.normalized_sensitivities),'FaceColor',[0.75,0.75,0.75]);
ax = gca;
ax.XTickLabel = SensitivityTable.CaseLabels;
set(ax,'XTickLabelRotation',90);

% caselabels = {'Base'; 'A'; 'B'; 'T'; 'E_r'; 'E_p'; 'NEPC 1';...
%     'NEPC 2'; 'NEPC 3'; 'NEPC 4'; 'NEPC 5';};
% bar(normalized_sensitivities,'FaceColor',[0.75,0.75,0.75]);
% set(gca,'xticklabel',caselabels);
% if ylimits == 1
%     ylim([lowlim,uplim]);
% end
title('Parameter Sensitivities (10% Change)');
% xlabel('Case');
ylabel('Normalized Sensitivity');
% ytickformat('percentage')
ylim([-0.05 1.03]);

hold on

er = errorbar(1:10,abs(SensitivityTable.normalized_sensitivities),SensitivityTable.normalized_stderror,SensitivityTable.normalized_stderror);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

hold off


figure(3);
bar(abs(SensitivityTable.normalized_sensitivities(6:end)),'FaceColor',[0.75,0.75,0.75]);
ax = gca;
ax.XTickLabel = SensitivityTable.CaseLabels(6:end);
set(ax,'XTickLabelRotation',90);
title('Parameter Sensitivities (10% Change)');
ylabel('Normalized Sensitivity');
hold on

er = errorbar(1:5,abs(SensitivityTable.normalized_sensitivities(6:end)),...
    SensitivityTable.normalized_stderror(6:end),SensitivityTable.normalized_stderror(6:end));
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

hold off


end