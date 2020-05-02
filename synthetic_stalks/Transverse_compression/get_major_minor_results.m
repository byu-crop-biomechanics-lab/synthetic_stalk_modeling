function [Results_new_x,Results_new_y,percents_x,percents_y] = get_major_minor_results(Data)
% FILENAME: get_sensitivity_resultsV2.m
% AUTHOR: Ryan Larson
% DATE: 2/19/2020
%
% PURPOSE: 
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

for j = 1:length(ResultsCell)
    Results = cell2mat(ResultsCell(j));
    rows = length(unique(Results(:,1)));
    cols = 2;
    Results_temp_x = NaN(rows,cols);
    Results_temp_y = NaN(rows,cols);

    for i = 1:size(Results,1)
        if i == 1
            row = 1;
        elseif Results(i,1) ~= Results(i-1,1)
            row = row + 1;
        end
    %     row = round(Results(i,1),0);
        col = round(Results(i,2),0) + 1;
        
        Results_temp_y(row,col) = abs(Results(i,4));
        if Results(i,2) == 1
            Results_temp_x(row,col) = abs(Results(i,3));
        else
            Results_temp_x(row,col) = abs(Results(i,4));
        end
        
%         if Results(i,2) == 0
%             Results_temp_x(row,col) = abs(Results(i,4));
%             Results_temp_x(row,col) = abs(Results(i,3));
%         else
%             Results_temp_y(row,col) = abs(Results(i,4));
%             Results_temp_y(row,col) = abs(Results(i,4));
%         end
    end

    if j == 1
        Results_new_x = Results_temp_x;
        Results_new_y = Results_temp_y;
    else
        Results_new_x = [Results_new_x; Results_temp_x];
        Results_new_y = [Results_new_y; Results_temp_y];
    end
end


percents_x = zeros(size(Results_new_x));
percents_y = zeros(size(Results_new_y));

for i = 1:size(Results_new_x,1)
    for j = 1:size(Results_new_x,2)
        percents_x(i,j) = (Results_new_x(i,j)/Results_new_x(i,1))*100;
    end
end

for i = 1:size(Results_new_y,1)
    for j = 1:size(Results_new_y,2)
        percents_y(i,j) = (Results_new_y(i,j)/Results_new_y(i,1))*100;
    end
end

XNans = isnan(percents_x(:,2));
YNans = isnan(percents_y(:,2));

percents_x(any(isnan(percents_x),2),:) = [];
percents_y(any(isnan(percents_y),2),:) = [];
percents_y(any(percents_y>200,2),:) = [];

error_x = percents_x - 100;
error_y = percents_y - 100;

figure(1);
boxplot(error_x(:,2));
title('Comparing x Component of Major Loading to Minor Loading');
ytickformat('percentage');
ylabel('Percent Error');

figure(2);
boxplot(error_y(:,2));
title('Comparing y Component of Major Loading to Minor Loading');
ytickformat('percentage');
ylabel('Percent Error');

end