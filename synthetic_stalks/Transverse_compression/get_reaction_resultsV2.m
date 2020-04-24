function [cumulative,cumulative_errbar,cumulative_err,individual,individual_errbar,individual_err,percents,percent_err,Results_new] = get_reaction_resultsV2(stalknums,numNEPCs,AllGoodReactionData)
% FILENAME: get_reaction_resultsV2.m
% AUTHOR: Ryan Larson
% DATE: 1/2020
%
% PURPOSE: Calculate the reaction results found in Results cell array (as
% of 4/8/2020, this is found in Results_AllPCA.mat).
% 
% 
% INPUTS:
%       AllGoodReactionData: 
%       
% OUTPUTS:
%       cumulative: 
% 
%       cumulative_errbar:
% 
%       cumulative_err:
% 
%       individual:
% 
%       individual_errbar:
% 
%       individual_err:
% 
%       percents: 
% 
%       percent_err:
% 
%       Results_new: 
%
%
% NOTES: 
%       
%       
% 
% 
% VERSION HISTORY:
% V1 - Designed to work with data from single slice locations only.
% V2 - Designed to work with ResultsCell format for combining results of
% multiple slice locations at once.
% V3 - 
%
% -------------------------------------------------------------------------

close;
set(0,'DefaultFigureWindowStyle','docked');
load(AllGoodReactionData,'ResultsCell');

% CONVERT FROM MICROMETER SCALE TO MILLIMETER SCALE FOR LOOKING AT ACTUAL
% VALUES

rows = length(stalknums);
cols = 1 + 2*numNEPCs;

ylimits = 0;
lowlim = 94;
uplim = 102;

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

%% "99% of the way there" method
% Convert to percentages of the original response
percents = zeros(size(Results_new));

for i = 1:size(Results_new,1)
    for j = 1:size(Results_new,2)
        percents(i,j) = (Results_new(i,j)/Results_new(i,1))*100;
    end
end

% Examine percentage results
avg = nanmean(percents);
stdev = nanstd(percents);
stderror = zeros(size(stdev));
n = sum(~isnan(percents),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev)
    stderror(i) = stdev(i)/sqrt(n(i));
end

%% Error method
% Create error population by shifting percentage values by 100%
percent_err = percents - 100;
percent_labels = strings(size(percent_err));

% Create labels according to the number of principal components used in
% the study (cumulative cases followed by remaining individual cases)
all_labels = strings(1,(1+2*numNEPCs));
all_labels = ["Real","Ellipse","Ellipse + PC 1"];
for i = 2:numNEPCs
    addlabel = "Ellipse + PCs 1-" + num2str(i);
    all_labels(1,i+2) = addlabel;
end

PCcase = 1;
for i = (3+numNEPCs):(1+2*numNEPCs)
    PCcase = PCcase + 1;
    addlabel = "Ellipse + PC " + num2str(PCcase);
    all_labels(1,i) = addlabel;
end

% % Fill percent_labels with copies of all_labels in the rows. This is the
% % form needed to do the boxplot labeling later on.
% percent_labels_test = percent_labels;
% for i = 1:size(percent_err,1)
%     percent_labels_test(i,:) = all_labels;
% end
% 
% all_labels_test = {'Real','Ellipse','Ellipse + PC 1','Ellipse + PCs 1-2',...
%     'Ellipse + PCs 1-3','Ellipse + PCs 1-4','Ellipse + PCs 1-5',...
%     'Ellipse + PC 2','Ellipse + PC 3','Ellipse + PC 4','Ellipse + PC 5'};
    
for i = 1:size(percent_err,1)
%     percent_labels(i,:) = {'Real','Ellipse','Ellipse + PC 1','Ellipse + PCs 1-2',...
%         'Ellipse + PCs 1-3','Ellipse + PCs 1-4','Ellipse + PCs 1-5',...
%         'Ellipse + PC 2','Ellipse + PC 3','Ellipse + PC 4','Ellipse + PC 5'};
    percent_labels(i,:) = all_labels;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examine percentage results (mostly used for checking outputs)
avg_err = nanmean(percent_err);
med_err = nanmedian(percent_err);
stdev_err = nanstd(percent_err);
stderror_err = zeros(size(stdev_err));
n_err = sum(~isnan(percent_err),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev_err)
    stderror_err(i) = stdev_err(i)/sqrt(n_err(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% Prepare data for boxplots
percent_box = [];
percent_boxlabels = [];
for i = 2:7
    percent_box = [percent_box; percent_err(:,i)];
    percent_boxlabels = [percent_boxlabels; percent_labels(:,i)];
end










% % % % Prepare data for boxplots
% % % percent_box = [];
% % % percent_boxlabels = [];
% % % for i = 2:(2+numNEPCs)
% % %     percent_box = [percent_box; percent_err(:,i)];
% % %     percent_boxlabels = [percent_boxlabels; percent_labels(:,i)];
% % % end

percent_box_ind = [];
percent_boxlabels_ind = [];
individual = linspace((3+numNEPCs),(1+2*numNEPCs),(numNEPCs-1));
ind_cases = [2 3 individual];
for i = ind_cases
    percent_box_ind = [percent_box_ind; percent_err(:,i)];
    percent_boxlabels_ind = [percent_boxlabels_ind; percent_labels(:,i)];
end


% % Calculate relative error for boxplots
% errdiffs = zeros(size(percent_err,1),10);
% errlabels = zeros(size(errdiffs));
% 
% for i = 1:size(errdiffs,1)
%     for j = 1:size(errdiffs,2)
%         errdiffs(i,j) = percent_err(i,j+1) - percent_err(i,j);
%         errlabels(i,j) = j;
%     end
% end
% 
% % Prepare for boxplots of the relative error
% boxerrdiffs = [];
% boxerrlabels = [];
% for i = 1:6
%     boxerrdiffs = [boxerrdiffs; errdiffs(:,i)];
%     boxerrlabels = [boxerrlabels; errlabels(:,i)];
% end




%% Plots (error method)
caselabels_cumulative = all_labels{1:(2+numNEPCs)};
caselabels_cumulative = caselabels_cumulative';

cumulative = avg(1:(2+numNEPCs));
cumulative_errbar = stderror(1:(2+numNEPCs));

caselabels_individual = {all_labels{1}, all_labels{2}, all_labels{3}, all_labels{(3+numNEPCs):end}};
caselabels_individual = caselabels_individual';
indices = [1 2 3 8 9 10 11];

individual = avg(indices);
individual_errbar = stderror(indices);

cumulative_err = avg_err(1:(2+numNEPCs));
% cumulative_errbar_err = stderror_err(1:7);

% hold on
% 
% er = errorbar(1:7,avg_err(1:7),stderror_err(1:7),stderror_err(1:7));
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% er.LineWidth = 0.5;
% 
% hold off

% % Individual NEPC bar chart (error method)
% caselabels_individual = {'Real'; 'Ellipse'; 'Ellipse + PC 1'; 'Ellipse + PC 2';...
%     'Ellipse + PC 3'; 'Ellipse + PC 4'; 'Ellipse + PC 5'};
% indices = [1 2 3 8 9 10 11];

individual_err = avg_err(indices);
% individual_errbar_err = stderror_err(indices);

% hold on
% 
% er = errorbar(1:7,avg_err(indices),stderror_err(indices),stderror_err(indices));
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% er.LineWidth = 0.5;


%% Boxplots (cumulative, percent error and relative percent error)
% Percent error
% Cumulative PC cases
figure(1);

boxplot(percent_box,percent_boxlabels,'Notch','on','symbol','');
ylim([-4,2]); % NEED TO SET YLIM ACCORDING TO WHERE OUTLIERS ARE (IGNORE THEM)
set(gca,'YTick',-4:0.5:2,'XTickLabelRotation',-30);
ytickformat('percentage');
ylabel('Error');
hold on
yline(0);
hold off

% Individual PC cases
figure(2);
boxplot(percent_box_ind,percent_boxlabels_ind,'Notch','on','symbol','');
ylim([-6.5,2.5]); % NEED TO SET YLIM ACCORDING TO WHERE OUTLIERS ARE (IGNORE THEM)
set(gca,'YTick',-6.5:0.5:2.5,'XTickLabelRotation',-30);
ytickformat('percentage');
ylabel('Error');
hold on
yline(0);
hold off

end