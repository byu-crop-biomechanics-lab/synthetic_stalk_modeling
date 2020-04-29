function [SensitivityTable,Stiffnesses] = get_sensitivity_resultsV2(stalknums,displacement,numNEPCs,pct_change,Data)
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
set(0,'DefaultFigureWindowStyle','docked');
load(Data,'ResultsCell');


% Define the size of the data arrays Results_new and Stiffnesses, which are
% organized with cases in the columns and slices in the rows
rows = length(stalknums);
cols = 1 + 2*numNEPCs;
% rows = 70;
% cols = 11;

Results_new = NaN(rows,cols);



nslices = length(ResultsCell);
for j = 1:nslices
    Results_temp = NaN(rows,cols);
    Results = cell2mat(ResultsCell(j));

    for i = 1:size(Results,1)
        if i == 1
            row = 1;
        elseif Results(i,1) ~= Results(i-1,1)
            row = row + 1;
        end
        
        col = round(Results(i,2),0) + 1;

        Results_temp(row,col) = abs(Results(i,4));
    end

    if j == 1
        Results_new = Results_temp;
    else
        Results_new = [Results_new; Results_temp];
    end
end

% Calculate the actual transverse stiffnesses in Pascals (displacement
% needs to be in meters, even though it was micrometers in FEA)
Stiffnesses = Results_new/displacement;


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

% Local sensitivity bar chart
caselabels = strings(1,(6+numNEPCs));
caselabels(1,1:6) = ["Base", "Major Diameter", "Minor Diameter",...
    "Rind Thickness", "Rind Modulus", "Pith Modulus"];
for i = 1:numNEPCs
    addlabel = "PC " + num2str(i);
    caselabels(1,i+6) = addlabel;
end
CaseLabels = categorical(caselabels);
CaseLabels = CaseLabels';

SensitivityTable = table(avg_sensitivities,CaseLabels,stderror,normalized_sensitivities,normalized_stderror,abs_avg_sensitivities);

% Sort table values by sensitivities
SensitivityTable = sortrows(SensitivityTable,'abs_avg_sensitivities','descend');
SensitivityTable(end,:) = [];

% Normalized sensitivity bar chart
figure(1);
bar(abs(SensitivityTable.normalized_sensitivities),'FaceColor',[0.75,0.75,0.75]);
ax = gca;
ax.XTickLabel = SensitivityTable.CaseLabels;
set(ax,'XTickLabelRotation',90);
ylabel('Normalized Sensitivity');
ylim([-0.05 1.05]);

hold on
er = errorbar(1:(5+numNEPCs),abs(SensitivityTable.normalized_sensitivities),SensitivityTable.normalized_stderror,SensitivityTable.normalized_stderror);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;
hold off

% Bar chart that shows just the smallest sensitivities (somewhat hard-coded
% after noting that the PC parameters all 
figure(2);
bar(abs(SensitivityTable.normalized_sensitivities(6:end)),'FaceColor',[0.75,0.75,0.75]);
ax = gca;
ax.XTickLabel = SensitivityTable.CaseLabels(6:end);
set(ax,'XTickLabelRotation',90);
ylabel('Normalized Sensitivity');

hold on
er = errorbar(1:numNEPCs,abs(SensitivityTable.normalized_sensitivities(6:end)),...
    SensitivityTable.normalized_stderror(6:end),SensitivityTable.normalized_stderror(6:end));
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;
hold off


end