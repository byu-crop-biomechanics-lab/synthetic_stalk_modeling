function get_material_sensitivity(Results,slices_unique,pct_change)
% FILENAME: get_sensitivity_resultsV2.m
% AUTHOR: Ryan Larson
% DATE: 2/3/2020
%
% PURPOSE: Wrap the majority of the data production process into a single
% script
% 
% 
% INPUTS:
%       Results: Saved in Material_Sensitivity_Data.mat
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

%% Prepare data into an array
nunique = size(slices_unique,1);

slices_unique = sortrows(slices_unique);

% Working array to hold values. Columns are cases, groups of 4 rows at a
% time corresponding to the groups
Results_new = NaN((nunique),6);

for i = 1:length(Results.Slice)
    % Get row in Results_new that corresponds to the data in Results
    A = slices_unique == [Results.Slice(i) Results.Stalk(i)];
    
    for j = 1:nunique
        if A(j,1) == 1 && A(j,2) == 1
            idx = j;
            break
        end        
    end
    
    row = 4*(idx-1) + Results.MatGroup(i);
    
    col = Results.Case(i) + 1;
    
    Results_new(row,col) = Results.RFy(i);
    
end

%% Calculate sensitivities
sensitivities = Results_new;

for i = 1:size(sensitivities,1)
    for j = 2:size(sensitivities,2)
        sensitivities(i,j) = get_sensitivity(sensitivities(i,1),sensitivities(i,j),pct_change);
    end
    sensitivities(i,1) = 0;
end

% Sort sensitivities into groups based on material sampling point
rows_minRminP = 1:4:size(sensitivities,1);
rows_minRmaxP = 2:4:size(sensitivities,1);
rows_maxRmaxP = 3:4:size(sensitivities,1);
rows_maxRminP = 4:4:size(sensitivities,1);

minRminP = [];
minRmaxP = [];
maxRmaxP = [];
maxRminP = [];

for row = rows_minRminP
    minRminP = [minRminP; sensitivities(row,:)];
end

for row = rows_minRmaxP
    minRmaxP = [minRmaxP; sensitivities(row,:)];
end

for row = rows_maxRmaxP
    maxRmaxP = [maxRmaxP; sensitivities(row,:)];
end

for row = rows_maxRminP
    maxRminP = [maxRminP; sensitivities(row,:)];
end

A = minRminP;
B = minRmaxP;
C = maxRmaxP;
D = maxRminP;

%% Paired t-tests
% Comparison 1: AB
diffAB = A-B;
meansAB = mean(diffAB);
stdAB = std(diffAB);
[H_AB0,P_AB0,CI_AB0,STATS_AB0] = ttest(A(:,1),B(:,1));
[H_AB1,P_AB1,CI_AB1,STATS_AB1] = ttest(A(:,2),B(:,2));
[H_AB2,P_AB2,CI_AB2,STATS_AB2] = ttest(A(:,3),B(:,3));
[H_AB3,P_AB3,CI_AB3,STATS_AB3] = ttest(A(:,4),B(:,4));
[H_AB4,P_AB4,CI_AB4,STATS_AB4] = ttest(A(:,5),B(:,5));
[H_AB5,P_AB5,CI_AB5,STATS_AB5] = ttest(A(:,6),B(:,6));


% Comparison 2: BC
diffBC = B-C;
meansBC = mean(diffBC);
stdBC = std(diffBC);
[H_BC0,P_BC0,CI_BC0,STATS_BC0] = ttest(B(:,1),C(:,1));
[H_BC1,P_BC1,CI_BC1,STATS_BC1] = ttest(B(:,2),C(:,2));
[H_BC2,P_BC2,CI_BC2,STATS_BC2] = ttest(B(:,3),C(:,3));
[H_BC3,P_BC3,CI_BC3,STATS_BC3] = ttest(B(:,4),C(:,4));
[H_BC4,P_BC4,CI_BC4,STATS_BC4] = ttest(B(:,5),C(:,5));
[H_BC5,P_BC5,CI_BC5,STATS_BC5] = ttest(B(:,6),C(:,6));

% Comparison 3: CD
diffCD = C-D;
meansCD = mean(diffCD);
stdCD = std(diffCD);
[H_CD0,P_CD0,CI_CD0,STATS_CD0] = ttest(C(:,1),D(:,1));
[H_CD1,P_CD1,CI_CD1,STATS_CD1] = ttest(C(:,2),D(:,2));
[H_CD2,P_CD2,CI_CD2,STATS_CD2] = ttest(C(:,3),D(:,3));
[H_CD3,P_CD3,CI_CD3,STATS_CD3] = ttest(C(:,4),D(:,4));
[H_CD4,P_CD4,CI_CD4,STATS_CD4] = ttest(C(:,5),D(:,5));
[H_CD5,P_CD5,CI_CD5,STATS_CD5] = ttest(C(:,6),D(:,6));

% Comparison 4: DA
diffDA = D-A;
meansDA = mean(diffDA);
stdDA = std(diffDA);
[H_DA0,P_DA0,CI_DA0,STATS_DA0] = ttest(D(:,1),A(:,1));
[H_DA1,P_DA1,CI_DA1,STATS_DA1] = ttest(D(:,2),A(:,2));
[H_DA2,P_DA2,CI_DA2,STATS_DA2] = ttest(D(:,3),A(:,3));
[H_DA3,P_DA3,CI_DA3,STATS_DA3] = ttest(D(:,4),A(:,4));
[H_DA4,P_DA4,CI_DA4,STATS_DA4] = ttest(D(:,5),A(:,5));
[H_DA5,P_DA5,CI_DA5,STATS_DA5] = ttest(D(:,6),A(:,6));

% Comparison 5: AC
diffAC = A-C;
meansAC = mean(diffAC);
stdAC = std(diffAC);
[H_AC0,P_AC0,CI_AC0,STATS_AC0] = ttest(A(:,1),C(:,1));
[H_AC1,P_AC1,CI_AC1,STATS_AC1] = ttest(A(:,2),C(:,2));
[H_AC2,P_AC2,CI_AC2,STATS_AC2] = ttest(A(:,3),C(:,3));
[H_AC3,P_AC3,CI_AC3,STATS_AC3] = ttest(A(:,4),C(:,4));
[H_AC4,P_AC4,CI_AC4,STATS_AC4] = ttest(A(:,5),C(:,5));
[H_AC5,P_AC5,CI_AC5,STATS_AC5] = ttest(A(:,6),C(:,6));

% Comparison 6: BD
diffBD = B-D;
meansBD = mean(diffBD);
stdBD = std(diffBD);
[H_BD0,P_BD0,CI_BD0,STATS_BD0] = ttest(B(:,1),D(:,1));
[H_BD1,P_BD1,CI_BD1,STATS_BD1] = ttest(B(:,2),D(:,2));
[H_BD2,P_BD2,CI_BD2,STATS_BD2] = ttest(B(:,3),D(:,3));
[H_BD3,P_BD3,CI_BD3,STATS_BD3] = ttest(B(:,4),D(:,4));
[H_BD4,P_BD4,CI_BD4,STATS_BD4] = ttest(B(:,5),D(:,5));
[H_BD5,P_BD5,CI_BD5,STATS_BD5] = ttest(B(:,6),D(:,6));


%% Export data
ABstats =  {H_AB0,P_AB0,CI_AB0,STATS_AB0;
            H_AB1,P_AB1,CI_AB1,STATS_AB1;
            H_AB2,P_AB2,CI_AB2,STATS_AB2;
            H_AB3,P_AB3,CI_AB3,STATS_AB3;
            H_AB4,P_AB4,CI_AB4,STATS_AB4;
            H_AB5,P_AB5,CI_AB5,STATS_AB5};
BCstats =  {H_BC0,P_BC0,CI_BC0,STATS_BC0;
            H_BC1,P_BC1,CI_BC1,STATS_BC1;
            H_BC2,P_BC2,CI_BC2,STATS_BC2;
            H_BC3,P_BC3,CI_BC3,STATS_BC3;
            H_BC4,P_BC4,CI_BC4,STATS_BC4;
            H_BC5,P_BC5,CI_BC5,STATS_BC5};
CDstats =  {H_CD0,P_CD0,CI_CD0,STATS_CD0;
            H_CD1,P_CD1,CI_CD1,STATS_CD1;
            H_CD2,P_CD2,CI_CD2,STATS_CD2;
            H_CD3,P_CD3,CI_CD3,STATS_CD3;
            H_CD4,P_CD4,CI_CD4,STATS_CD4;
            H_CD5,P_CD5,CI_CD5,STATS_CD5};
DAstats =  {H_DA0,P_DA0,CI_DA0,STATS_DA0;
            H_DA1,P_DA1,CI_DA1,STATS_DA1;
            H_DA2,P_DA2,CI_DA2,STATS_DA2;
            H_DA3,P_DA3,CI_DA3,STATS_DA3;
            H_DA4,P_DA4,CI_DA4,STATS_DA4;
            H_DA5,P_DA5,CI_DA5,STATS_DA5};
ACstats =  {H_AC0,P_AC0,CI_AC0,STATS_AC0;
            H_AC1,P_AC1,CI_AC1,STATS_AC1;
            H_AC2,P_AC2,CI_AC2,STATS_AC2;
            H_AC3,P_AC3,CI_AC3,STATS_AC3;
            H_AC4,P_AC4,CI_AC4,STATS_AC4;
            H_AC5,P_AC5,CI_AC5,STATS_AC5};
BDstats =  {H_BD0,P_BD0,CI_BD0,STATS_BD0;
            H_BD1,P_BD1,CI_BD1,STATS_BD1;
            H_BD2,P_BD2,CI_BD2,STATS_BD2;
            H_BD3,P_BD3,CI_BD3,STATS_BD3;
            H_BD4,P_BD4,CI_BD4,STATS_BD4;
            H_BD5,P_BD5,CI_BD5,STATS_BD5};

FolderName = pwd;
FileName = 'MaterialSensitivityResults.mat';
SaveFile = fullfile(FolderName, FileName);
save(SaveFile,'slices_unique','A','B','C','D','diffAB','meansAB','stdAB',...
    'diffBC','meansBC','stdBC','diffCD','meansCD','stdCD','diffDA',...
    'meansDA','stdDA','diffAC','meansAC','stdAC','diffBD','meansBD',...
    'stdBD','ABstats','BCstats','CDstats','DAstats','ACstats','BDstats');


end