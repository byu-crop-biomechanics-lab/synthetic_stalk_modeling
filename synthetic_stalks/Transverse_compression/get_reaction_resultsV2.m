function [cumulative,cumulative_errbar,cumulative_err,individual,individual_errbar,individual_err,percents,percent_err] = get_reaction_resultsV2(AllGoodReactionData)
% get_reaction_resultsV2.m: Take the text file structure of the NEPC
% finite-element results and reorganize the data into an array, with the
% rows being the cross-section number and the columns being the case
% number. Cases 0-6 are used.
% clear; close;

close;
load(AllGoodReactionData,'ResultsCell');

% CONVERT FROM MICROMETER SCALE TO MILLIMETER SCALE FOR LOOKING AT ACTUAL
% VALUES

rows = 50;
cols = 11;

ylimits = 0;
lowlim = 94;
uplim = 102;

Results_new = NaN(rows,cols);

for j = 1:4
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
percent_labels = cell(size(percent_err));

for i = 1:size(percent_err,1)
    percent_labels(i,:) = {'Real','Ellipse','NEPC 1','NEPC 1-2','NEPC 1-3','NEPC 1-4','NEPC 1-5','NEPC 2','NEPC 3','NEPC 4','NEPC 5'};
end

% Examine percentage results
avg_err = nanmean(percent_err);
stdev_err = nanstd(percent_err);
stderror_err = zeros(size(stdev_err));
n_err = sum(~isnan(percent_err),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev_err)
    stderror_err(i) = stdev_err(i)/sqrt(n_err(i));
end

% Prepare data for boxplots
percent_box = [];
percent_boxlabels = [];
for i = 2:7
    percent_box = [percent_box; percent_err(:,i)];
    percent_boxlabels = [percent_boxlabels; percent_labels(:,i)];
end


% Calculate relative error for boxplots
errdiffs = zeros(size(percent_err,1),10);
errlabels = zeros(size(errdiffs));

for i = 1:size(errdiffs,1)
    for j = 1:size(errdiffs,2)
        errdiffs(i,j) = percent_err(i,j+1) - percent_err(i,j);
        errlabels(i,j) = j;
    end
end

% Prepare for boxplots of the relative error
boxerrdiffs = [];
boxerrlabels = [];
for i = 1:6
    boxerrdiffs = [boxerrdiffs; errdiffs(:,i)];
    boxerrlabels = [boxerrlabels; errlabels(:,i)];
end




%% Plots
% Cumulative case bar chart ("99% of the way there" method)
% BREAK HERE IF YOU WANT TO GET INTERMEDIATE PERCENT VALUES
figure(1);
caselabels_cumulative = {'Real'; 'Ellipse'; 'NEPC 1'; 'NEPC 1-2'; 'NEPC 1-3'; 'NEPC 1-4'; 'NEPC 1-5'};
bar(avg(1:7),'FaceColor',[0.75,0.75,0.75]);
set(gca,'xticklabel',caselabels_cumulative);
if ylimits == 1
    ylim([lowlim,uplim]);
end
title('NEPC Response Progression (Cumulative)');
xlabel('Case');
ylabel('Reaction Force (% of Real Response)');

cumulative = avg(1:7);
cumulative_errbar = stderror(1:7);

hold on

er = errorbar(1:7,avg(1:7),stderror(1:7),stderror(1:7));
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

% % Define noisy region as +/- 0.2% of the real cross-section response
% yline(99.8,':k');
% yline(100.2,':k');

hold off

% Individual NEPC bar chart ("99% of the way there" method)
caselabels_individual = {'Real'; 'Ellipse'; 'NEPC 1'; 'NEPC 2'; 'NEPC 3'; 'NEPC 4'; 'NEPC 5'};
indices = [1 2 3 8 9 10 11];
figure(2);
bar(avg(indices),'FaceColor',[0.75,0.75,0.75]);
set(gca,'xticklabel',caselabels_individual);
if ylimits == 1
    ylim([lowlim,uplim]);
end
title('NEPC Response Progression (Individual)');
xlabel('Case');
ylabel('Reaction Force (% of Real Response)');

individual = avg(indices);
individual_errbar = stderror(indices);

hold on

er = errorbar(1:7,avg(indices),stderror(indices),stderror(indices));
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

% % Define noisy region as +/- 0.2% of the real cross-section response
% yline(99.8,':k');
% yline(100.2,':k');

hold off

% ============================
%         Error method
% ============================

% Cumulative case bar chart (error method)
figure(3);
caselabels_cumulative = {'Real'; 'Ellipse'; 'NEPC 1'; 'NEPC 1-2'; 'NEPC 1-3'; 'NEPC 1-4'; 'NEPC 1-5'};
bar(avg_err(1:7),'FaceColor',[0.75,0.75,0.75]);
set(gca,'xticklabel',caselabels_cumulative);
if ylimits == 1
    ylim([lowlim,uplim]);
end
title('NEPC Response Progression (Cumulative)');
xlabel('Case');
ylabel('Reaction Force Error (% of Real Response)');

cumulative_err = avg_err(1:7);
cumulative_errbar_err = stderror_err(1:7);

hold on

er = errorbar(1:7,avg_err(1:7),stderror_err(1:7),stderror_err(1:7));
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;

% % Define noisy region as +/- 0.2% of the real cross-section response
% yline(99.8,':k');
% yline(100.2,':k');

hold off

% Individual NEPC bar chart (error method)
caselabels_individual = {'Real'; 'Ellipse'; 'NEPC 1'; 'NEPC 2'; 'NEPC 3'; 'NEPC 4'; 'NEPC 5'};
indices = [1 2 3 8 9 10 11];
figure(4);
bar(avg_err(indices),'FaceColor',[0.75,0.75,0.75]);
set(gca,'xticklabel',caselabels_individual);
if ylimits == 1
    ylim([lowlim,uplim]);
end
title('NEPC Response Progression (Individual)');
xlabel('Case');
ylabel('Reaction Force Error (% of Real Response)');

individual_err = avg_err(indices);
individual_errbar_err = stderror_err(indices);

hold on

er = errorbar(1:7,avg_err(indices),stderror_err(indices),stderror_err(indices));
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.5;


%% Boxplots (cumulative, percent error and relative percent error)
% Percent error
figure(5);
boxplot(percent_box,percent_boxlabels,'Notch','on','symbol','');
ylim([-4.5,2]);
set(gca,'YTick',-4.5:0.5:2);
ytickformat('percentage');
title('Reaction Error Distributions');
ylabel('Error');
hold on
yline(0);
hold off

% Relative percent error
figure(6);
% boxplot(boxerrdiffs,boxerrlabels,'Notch','on');
boxplot(boxerrdiffs,boxerrlabels,'Notch','on','symbol','');
ylim([-4.5,3.5]);
set(gca,'YTick',-4.5:0.5:3.5);
ytickformat('percentage');
title('Reaction Error Distributions (Relative)');
ylabel('Relative Error');
hold on
yline(0);
hold off

end