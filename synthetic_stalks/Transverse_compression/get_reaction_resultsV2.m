function [cumulative_err,individual_err,percent_err,Stiffnesses] = get_reaction_resultsV2(stalknums,displacement,numNEPCs,AllGoodReactionData)
% FILENAME: get_reaction_resultsV2.m
% AUTHOR: Ryan Larson
% DATE: 1/2020
%
% PURPOSE: Calculate the reaction results found in Results cell array (as
% of 4/8/2020, this is found in Results_AllPCA.mat).
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
%       AllGoodReactionData: Input .mat file with the reaction results
%       data, organized in a Results cell array (as of 4/8/2020, this is
%       found in Results_AllPCA.mat).
%       
% OUTPUTS: 
%       cumulative_err: The median percent errors for the cumulative cases.
% 
%       individual_err: The median percent errors for the individual cases.
% 
%       percent_err: An array of percent error values, where the columns
%       correspond to the model cases (i.e. Real, Ellipse, Ellipse + PC 1,
%       etc.) and the rows correspond to individual stalks.
% 
%       Stiffnesses: The calculated stiffness values, in Pascals
%
%
% NOTES: 
%       - All outputs can be calculated from the raw data in Results_new.
%       The median values for the cumulative cases and individual cases, as
%       well as the general data in percent_err, are provided for
%       convenience.
% -------------------------------------------------------------------------
% SUBROUTINES:
%   RoundToPoint5.m: Take an input double and round it to the nearest half.
% 
%   getBoxLims.m: Calculate roughly where the box plot whiskers will end so
%   appropriate ylim values can be chosen to show just the area around the
%   box plot when the outliers are ignored.
% 
% PSEUDO-CODE:
%   
% 
% -------------------------------------------------------------------------
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

% Define the size of the data arrays Results_new and Stiffnesses, which are
% organized with cases in the columns and slices in the rows
rows = length(stalknums);
cols = 1 + 2*numNEPCs;

% Create row vectors of the indices for the relevant cases in each
% experimental design (cumulative and individual). Ignore the first index,
% which is the real case, and will analytically have 0 percent error
% (however, the stiffness values might be of interest).
cumulative_indices = [2 linspace(3,2+numNEPCs,numNEPCs)];
individual_indices = [2 3 linspace((3+numNEPCs),(1+2*numNEPCs),(numNEPCs-1))];

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


%% "99% of the way there" method
% Note: All the calculations downstream rely on Results_new, which is the
% raw reaction force data. Since the applied displacement was linear, the
% percent error values will be the same whether we're operating on the
% reaction force data or the stiffness data.

% Convert to percentages of the original response
percents = zeros(size(Results_new));

for i = 1:size(Results_new,1)
    for j = 1:size(Results_new,2)
        percents(i,j) = (Results_new(i,j)/Results_new(i,1))*100;
    end
end

% Examine percentage results
stdev = nanstd(percents);
stderror = zeros(size(stdev));
n = sum(~isnan(percents),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev)
    stderror(i) = stdev(i)/sqrt(n(i));
end

%% Error method (what we're really interested in)
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
    
for i = 1:size(percent_err,1)
    percent_labels(i,:) = all_labels;
end

% Examine percentage results (mostly used for checking outputs)
med_err = nanmedian(percent_err);
stdev_err = nanstd(percent_err);
stderror_err = zeros(size(stdev_err));
n_err = sum(~isnan(percent_err),1); % row vector of total non-NaN elements in each column
for i = 1:length(stdev_err)
    stderror_err(i) = stdev_err(i)/sqrt(n_err(i));
end

% Prepare data for boxplots
percent_box = [];
percent_boxlabels = [];
for i = 2:(2+numNEPCs)
    percent_box = [percent_box; percent_err(:,i)];
    percent_boxlabels = [percent_boxlabels; percent_labels(:,i)];
end

percent_box_ind = [];
percent_boxlabels_ind = [];
individual = linspace((3+numNEPCs),(1+2*numNEPCs),(numNEPCs-1));
ind_cases = [2 3 individual];
for i = ind_cases
    percent_box_ind = [percent_box_ind; percent_err(:,i)];
    percent_boxlabels_ind = [percent_boxlabels_ind; percent_labels(:,i)];
end

%% Calculate outputs
% Median errors for each case in each experimental design
cumulative_err = med_err(cumulative_indices);
individual_err = med_err(individual_indices);


%% Box plots
% Get the uplim and lolim values for the cumulative cases (prevent the box
% plots from extending beyond the plot bounds)
uplimrow = zeros(1,1+numNEPCs);
lolimrow = zeros(1,1+numNEPCs);
% cumulative_indices = [2 linspace(3,2+numNEPCs,numNEPCs)];
for i = 1:length(uplimrow)
    [uplimrow(i),lolimrow(i)] = getBoxLims(percent_err(:,cumulative_indices(i)));
end

% Add a buffer between the calculated outer reach of the whiskers and the
% edge of the plot
uplim = RoundToPoint5(max(uplimrow) + 0.5);
lolim = RoundToPoint5(min(lolimrow) - 0.5);

% Box plot for cumulative case models
figure(1);
boxplot(percent_box,percent_boxlabels,'Notch','on','symbol','');
ylim([lolim,uplim]);
set(gca,'YTick',lolim:0.5:uplim,'XTickLabelRotation',-30);
ytickformat('percentage');
ylabel('Error');
hold on
yline(0);
hold off


% Get the uplim and lolim values for the individual cases (prevent the box
% plots from extending beyond the plot bounds)
uplimrow = zeros(1,1+numNEPCs);
lolimrow = zeros(1,1+numNEPCs);
% individual_indices = [2 3 linspace((3+numNEPCs),(1+2*numNEPCs),(numNEPCs-1))];
for i = 1:length(uplimrow)
    [uplimrow(i),lolimrow(i)] = getBoxLims(percent_err(:,individual_indices(i)));
end

% Add a buffer between the calculated outer reach of the whiskers and the
% edge of the plot
uplim = RoundToPoint5(max(uplimrow) + 0.5);
lolim = RoundToPoint5(min(lolimrow) - 0.5);

% Box plot for individual case models
figure(2);
boxplot(percent_box_ind,percent_boxlabels_ind,'Notch','on','symbol','');
ylim([lolim,uplim]); % NEED TO SET YLIM ACCORDING TO WHERE OUTLIERS ARE (IGNORE THEM)
set(gca,'YTick',lolim:0.5:uplim,'XTickLabelRotation',-30);
ytickformat('percentage');
ylabel('Error');
hold on
yline(0);
hold off

end

%% Localized functions
function [roundval] = RoundToPoint5(input)
% Round a number to a the closest 0.5

% Determine whether to round up or down, and then perform that operation
modval = mod(input,0.5);

if modval >= 0.25
    roundval = ceil(input*2)/2;
else
    roundval = floor(input*2)/2;
end

end


function [uplim,lolim] = getBoxLims(datacol)
% Take in a column vector of data going into a box plot and determine
% reasonable upper and lower limits for the whiskers. These won't be
% exactly what Matlab puts out, but it will be close.

Q1 = quantile(datacol,0.25);
Q3 = quantile(datacol,0.75);
IQR = Q3 - Q1;

uplim = Q3 + 1.5*IQR;
lolim = Q1 - 1.5*IQR;

end