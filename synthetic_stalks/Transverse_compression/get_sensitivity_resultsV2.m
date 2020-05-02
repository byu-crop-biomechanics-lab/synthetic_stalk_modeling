function [SensitivityTable,Stiffnesses] = get_sensitivity_resultsV2(stalknums,displacement,numNEPCs,pct_change,Data)
% FILENAME: get_sensitivity_resultsV2.m
% AUTHOR: Ryan Larson
% DATE: 2/3/2020
%
% PURPOSE: Calculate sensitivities and plot them for each variable
% 
% 
% INPUTS:
%       stalknums: The random stalk numbers that were chosen for model
%       production and subsequent FEA.
% 
%       displacement: The deflection that was applied to each model, IN
%       METERS.
% 
%       numNEPCs: The number of principal components used in the model
%       study.
% 
%       pct_change: The percentage to change each parameter value by
%       for the sensitivity study. Enter as a decimal (i.e. 10% would be
%       0.1). This value must be the same value used as an input to
%       TransverseSensitivityV1 for the data being analyzed.
% 
%       Data: A .mat file that contains a ResultsCell array, which is the
%       output of a Results.txt file organized into a single-row cell
%       array. As of 4/30/2020, this data file is
%       'Results_AllSensitivity.mat'.
% 
% OUTPUTS:
%       SensitivityTable: A table of sensitivity values for each parameter,
%       where the row numbers correspond to the case numbers (the zeroth
%       case is the base case, and therefore doesn't have a sensitivity).
%
%       Stiffnesses: An array of stiffness values (should be in Pascals).
%       The columns correspond to cases (including the zeroth case) and the
%       rows correspond to unique cross-sections.
% 
% NOTES: 
%     
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Load results data and create a blank array to reorganize data.
%   Iterate through all data and populate the array with all data.
%   With all data input, calculate stiffnesses for each cross section.
%   Calculate sensitivities for each parameter and populate an array.
%   Create bar chart for sensitivities
%       Sort table values by sensitivities.
%       Examine percentage data for box plotting (find median errors, find 
%       upper and lower limits, and add a buffer from edge of whiskers to
%       edge of plot).
%       Create the boxplot figures using populated labels.
% 
% 
% -------------------------------------------------------------------------
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

% Create Results_new with NaN elements so missing cases can be ignored in
% statistical calculations
Results_new = NaN(rows,cols);

% Iterate through the results for each slice location and reorganize the
% data to fit the structure of Results_new, with each case having a
% specific column and each row representing a unique slice.
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

%% Calculate sensitivities
sensitivities = NaN(size(Results_new));

% Each sensitivity is calculated relative to the results of the base case.
% The sensitivity of the column corresponding to the base case will be
% zero.
for i = 1:size(sensitivities,1)
    for j = 1:size(sensitivities,2)
        sensitivities(i,j) = ((Results_new(i,j) - Results_new(i,1))/Results_new(i,1));
    end
end

% Calculate average sensitivities and relevant statistics
avg_sensitivities = nanmean(sensitivities);
stdev = nanstd(sensitivities);  % Standard deviation
stderror = zeros(size(stdev));  
n = sum(~isnan(sensitivities),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev)
    stderror(i) = stdev(i)/sqrt(n(i));  % Standard error
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