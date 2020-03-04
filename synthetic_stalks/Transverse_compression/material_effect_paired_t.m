function material_effect_paired_t(Max_Min_Results)

load(Max_Min_Results);

%% Reorganize data
rows = 70;
cols = 11;

Results_new_Max = NaN(rows,cols);
Results_new_Min = NaN(rows,cols);

for j = 1:3
    Results_Max_temp = NaN(rows,cols);
    Results = cell2mat(Max(j));
    
    for i = 1:size(Results,1)
        if i == 1
            row = 1;
        elseif Results(i,1) ~= Results(i-1,1)
            row = row + 1;
        end
        
        col = round(Results(i,2),0) + 1;
    
        Results_Max_temp(row,col) = abs(Results(i,4));
    end
    
    
    Results_Min_temp = NaN(rows,cols);
    Results = cell2mat(Min(j));
    
    for i = 1:size(Results,1)
        if i == 1
            row = 1;
        elseif Results(i,1) ~= Results(i-1,1)
            row = row + 1;
        end
        
        col = round(Results(i,2),0) + 1;
    
        Results_Min_temp(row,col) = abs(Results(i,4));
    end
    
    if j == 1
        Results_new_Max = Results_Max_temp;
        Results_new_Min = Results_Min_temp;
    else
        Results_new_Max = [Results_new_Max; Results_Max_temp];
        Results_new_Min = [Results_new_Min; Results_Min_temp];
    end
    
    
end
    
    
    
    
%% Calculate data as percentages
Results_Max_pct = Results_new_Max;
Results_Min_pct = Results_new_Min;

for i = 1:size(Results_Max_pct,1)
    for j = 2:size(Results_Max_pct,2)
        Results_Max_pct(i,j) = 100*Results_Max_pct(i,j)/Results_Max_pct(i,1);        
    end
    Results_Max_pct(i,1) = 100;
end


for i = 1:size(Results_Min_pct,1)
    for j = 2:size(Results_Min_pct,2)
        Results_Min_pct(i,j) = 100*Results_Min_pct(i,j)/Results_Min_pct(i,1);        
    end
    Results_Min_pct(i,1) = 100;
end


MaxMinDiffs = Results_Max_pct - Results_Min_pct;

%% SAVE THE INTERMEDIATE DATA ARRAYS!!!!!!!


%% Get statistics of interest

meanpct = nanmean(MaxMinDiffs);

stdpct = nanstd(MaxMinDiffs);

% Individual case paired t-tests
[~,p0,~,~] = ttest(MaxMinDiffs(:,1))
[~,p1,~,~] = ttest(MaxMinDiffs(:,2))
[~,p2,~,~] = ttest(MaxMinDiffs(:,3))
[~,p3,~,~] = ttest(MaxMinDiffs(:,4))
[~,p4,~,~] = ttest(MaxMinDiffs(:,5))
[~,p5,~,~] = ttest(MaxMinDiffs(:,6))
[~,p6,~,~] = ttest(MaxMinDiffs(:,7))
[~,p7,~,~] = ttest(MaxMinDiffs(:,8))
[~,p8,~,~] = ttest(MaxMinDiffs(:,9))
[~,p9,~,~] = ttest(MaxMinDiffs(:,10))
[~,p10,~,~] = ttest(MaxMinDiffs(:,11))


pvals_all = [p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
pvals_detail = [p1,p2,p3,p4,p5,p6];

figure(1);
allcaselabels = {'Real'; 'Ellipse'; 'Ellipse + PC 1'; 'Ellipse + PCs 1-2';...
    'Ellipse + PCs 1-3'; 'Ellipse + PCs 1-4'; 'Ellipse + PCs 1-5';...
    'Ellipse + PC 2'; 'Ellipse + PC 3'; 'Ellipse + PC 4'; 'Ellipse + PC 5'};
bar(pvals_all,'FaceColor',[0.75,0.75,0.75]);
set(gca,'xticklabel',allcaselabels,'XTickLabelRotation',90);
xlabel('Geometric Case');
ylabel('P-value');
title('Paired T-test for Material Effects');


% All data combined paired t-test
AllDiffs = [];
for i = 2:7
    AllDiffs = [AllDiffs; MaxMinDiffs(:,i)];
end

[~,pAll,~,~] = ttest(AllDiffs)










end
