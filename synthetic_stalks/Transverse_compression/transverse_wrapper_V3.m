function transverse_wrapper_V3(stalknums,slicedist,material_method,group)
% FILENAME: transverse_wrapper.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: Wrap the majority of the data production process into a single
% script
% 
% 
% INPUTS:
%       stalknums - A vector of unique integers from 1 to 980 that determines
%       which stalks to sample from (use randperm(980,K) to choose K
%       unique integers from 1 to 980)
%       slicedist - A number value (can be integer or decimal) that
%       signifies the distance in millimeters from the nearest node on the
%       stalk at the chosen slice.
%       
% OUTPUTS:
%       - Several .mat files with variables saved from the steps in the
%       process. These are made available for troubleshooting purposes.
%       - Lots of Python scripts, corresponding to the stalks used and the
%       cases examined. The names do not contain the slice distance info.
%       11 scripts are created per stalk.
%
%
% NOTES: 
%       - The closest slice to slicedist for a given stalk will be chosen,
%       so when verifying results make sure to look at the distribution of
%       chosen slice positions when choosing extremely large magnitude
%       numbers.
%       - Don't go above a magnitude of 90 at all, since there are
%       very few stalks that have slices that far from the node. A value of
%       0 for slicedist will give the closest slice to the node.
%       - The PCA functions look at the first five NEPCs, so the range
%       should always cover at least 5 stalks (realistically, the range
%       should cover much more)
% 
% 
% VERSION HISTORY:
% V1 - Choose a continuous range of stalk numbers to sample
% V2 - Input a vector of unique random integers to determine which stalks
% to sample
% V3 - Added group input to number the output scripts for putting in a
% large folder of scripts for different purposes
%
% -------------------------------------------------------------------------

%% Process
% Load in the 360-point DCR version of Jared's table
load StalksDCR_360pts.mat
hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

dist_int = num2str(floor(abs(slicedist)));
deci = abs(slicedist) - floor(abs(slicedist));
dist_deci = num2str(deci);
dist_deci = erase(dist_deci,'0.');

if slicedist > 0
    if strcmp(dist_deci,'0')
        slicepos = strcat('Random_Above_',dist_int);
    else
        slicepos = strcat('Random_Above_',dist_int,'_',dist_deci);
    end
elseif slicedist < 0
    if strcmp(dist_deci,'0')
        slicepos = strcat('Random_Below_',dist_int);
    else
        slicepos = strcat('Random_Below_',dist_int,'_',dist_deci);
    end
else
    slicepos = strcat('Random_At_Node');
end

output_prefix = strcat(slicepos);

StalkNumsName = strcat(output_prefix,'_Stalks.mat');
FolderName = pwd;
SaveFile = fullfile(FolderName, StalkNumsName);
save(SaveFile,'stalknums');


% Gather all the cross-sections at the chosen slice distance
AllSectionsName = strcat(output_prefix,'_All980.mat');
ChooseSections('samedist',linspace(1,980,980),slicedist,Stalk_TableDCR,error_indices,npoints,AllSectionsName)
load(AllSectionsName);


%% Check to see if there is already a flip vector for the chosen distance
% Get file name to look for
searchslice = sprintf('%d',abs(slicedist));
if slicedist > 0
    FlippedOutputName = strcat('*Above_',searchslice,'_FLIPPED.mat');
    fstruct = dir(FlippedOutputName);
elseif slicedist < 0
    FlippedOutputName = strcat('*Below_',searchslice,'_FLIPPED.mat');
    fstruct = dir(FlippedOutputName);
else
    FlippedOutputName = strcat('*At_Node','_FLIPPED.mat');
    fstruct = dir(FlippedOutputName);        
end

if isempty(fstruct)
    FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
else
    FlippedOutputName = fstruct(1).name;
end

FlipName = strcat(output_prefix,'_flip_sections.mat');
% FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
if ~isfile(FlippedOutputName)
    disp('No flip index vector exists in the current folder. Create one now.');
    % Manually find the cross-sections that need to be flipped 180 degrees
%     FlipName = strcat(output_prefix,'_flip_sections.mat');
    find_flip_notches(AllSectionsName,FlipName)

    while 1
        fixes_needed = input('Does the flip vector need manual correction? Y/N ','s');
        switch fixes_needed
            case 'Y'
                load(FlipName);
                openvar('flip_sections');
                disp('Giving control to keyboard for manual editing of flip variable.');
                disp('Use dbcont command to exit keyboard editing mode.');
                keyboard;
                break
            case 'N'
                break

            otherwise
                disp('Not a recognized answer. Please try again.');
        end
    end


    % Flip the cross-sections that need to be flipped, according to the vector
    % of flip indicators
    % load(FlipName);
    FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
    flip_notches(FlipName,AllSectionsName,FlippedOutputName);
    
else
%     FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
    disp('A flip index vector for the chosen data has been found.');
end


load(FlippedOutputName);

% Make smaller chosen sections set
ChooseSectionsName = strcat(output_prefix,'_Sampled.mat');
ChooseSections('samedist',stalknums,slicedist,flippedTable,error_indices,npoints,ChooseSectionsName);

%% Check to see if there is already an ellipse fit for the chosen distance
% Get file name to look for
searchslice = sprintf('%d',slicedist);
if slicedist > 0
    searchname = strcat('*Above_',searchslice,'_AllGoodEllipses.mat');
    fstruct = dir(searchname);
    if isempty(fstruct)
        searchname = strcat('*Above_',searchslice,'_AllEllipses.mat');
    end
        
elseif slicedist < 0
    searchname = strcat('*Below_',searchslice,'_AllGoodEllipses.mat');
    fstruct = dir(searchname);
    if isempty(fstruct)
        searchname = strcat('*Below_',searchslice,'_AllEllipses.mat');
    end
    
else
    searchname = strcat('*At_Node','_AllGoodEllipses.mat');
    fstruct = dir(searchname);
    if isempty(fstruct)
        searchname = strcat('*At_Node','_AllEllipses.mat');
    end
    
end
    
fstruct = dir(searchname);
if isempty(fstruct)
    searchname = strcat(output_prefix,'_AllEllipses.mat');
else
    searchname = fstruct(1).name;
end



% If the ellipse data already exists, use the existing data instead of
% going through the process of tagging any bad fits
problem_indices = [];
chosen_problem_indices = [];
    
if ~isfile(searchname)
    disp('No ellipse data for this slice distance exists in the current folder. Create data now.');
    
    % Calculate the ellipse fits for all cross-sections at chosen distance
    AllEllipseName = strcat(output_prefix,'_AllEllipses.mat');
    ellipse_fitting_V2(FlippedOutputName,AllEllipseName);

    % Calculate the ellipse fits for the chosen set of cross-sections
    ChosenEllipseName = strcat(output_prefix,'_ChosenEllipses.mat');
    ellipse_fitting_V2(ChooseSectionsName,ChosenEllipseName);

    % Plot the ellipse fits and see if any of them have problems
    load(AllEllipseName);
    load(FlippedOutputName);

    for i = 1:size(flippedTable,1)
        i
        polarplot(ELLIPSE_T(i,:),ELLIPSE_R_ext(i,:),'LineWidth',2);
        hold on
        polarplot(ELLIPSE_T(i,:),ELLIPSE_R_int(i,:),'LineWidth',2);
        polarplot(ext_T(:,i),ext_Rho(:,i));
        polarplot(int_T(:,i),int_Rho(:,i));
        hold off
        s = input('Enter 1 if the ellipse fit has a problem: ');
        s
        if s == 1
            problem_indices = [problem_indices, i];
%             if any(ismember(i,stalknums))
%                 chosen_problem_indices = [chosen_problem_indices, i];
%             end
        else
            continue
        end    
    end    
    
    while 1
        fixes_needed = input('Does the ellipse problems vector need manual correction? Y/N ','s');
        switch fixes_needed
            case 'Y'
%                 load(FlipName);
                openvar('problem_indices');
                disp('Giving control to keyboard for manual editing of flip variable.');
                disp('Use dbcont command to exit keyboard editing mode.');
                keyboard;
                break
            case 'N'
                break

            otherwise
                disp('Not a recognized answer. Please try again.');
        end
    end
    
    

    % Determine if any of the problem cases are in the selected stalks
    problocs = ismember(problem_indices,stalknums);
    chosen_problem_indices = problem_indices(problocs);
    
    % Save problem_indices for later use, after correcting the indices
    ProblemEllipses = strcat(output_prefix,'_ProblemEllipses.mat');
    save(ProblemEllipses,'problem_indices','problocs','chosen_problem_indices');
    
    problem_indices
    chosen_problem_indices
    
    
else
    disp('Ellipse data for the chosen slice distance has been found.');
    AllEllipseName = strcat(output_prefix,'_AllEllipses.mat');
    
    ProblemEllipses = strcat(output_prefix,'_ProblemEllipses.mat');
    load(ProblemEllipses,'problem_indices','chosen_problem_indices');
    
    % Calculate the ellipse fits for the chosen set of cross-sections
    ChosenEllipseName = strcat(output_prefix,'_ChosenEllipses.mat');
    ellipse_fitting_V2(ChooseSectionsName,ChosenEllipseName);
end



NEPCName = strcat(output_prefix,'_PCA.mat');
MaterialsName = strcat(output_prefix, '_Materials.mat');

if isempty(problem_indices)
    % Run PCA on all cases
    PCA_ellipse_fits(AllEllipseName,NEPCName);
    
    % Create the Abaqus Python scripts
    create_cases(NEPCName,ChosenEllipseName,ChooseSectionsName,problem_indices,5,material_method,group,MaterialsName);
else
    % Remove the problem ellipses and then run PCA
    AllGoodEllipseFits = strcat(output_prefix,'_AllGoodEllipses.mat');
    ChosenGoodEllipseFits = strcat(output_prefix,'_ChosenGoodEllipses.mat');
    GoodStalkNumsName = strcat(output_prefix,'_GoodStalks.mat');

    % Remove the problem indices from the complete set of ellipse fits
    remove_problem_ellipses(AllEllipseName,stalknums,problem_indices,AllGoodEllipseFits,GoodStalkNumsName,'all')
    
    % Convert the chosen_problem_indices to values from 1 to 50 (or the
    % number of chosen stalks) and remove the corresponding bad sections
    % all at once
    translated_prob_indices = [];
    for i = 1:length(stalknums)
        if any(ismember(stalknums(i),chosen_problem_indices))
            translated_prob_indices = [translated_prob_indices,i];
        end
    end
    translated_prob_indices
    
    % Remove the problem indices that are within the chosen stalknums and
    % produce a new stalknums variable
    if ~isempty(translated_prob_indices)
        remove_problem_ellipses(ChosenEllipseName,stalknums,translated_prob_indices,ChosenGoodEllipseFits,GoodStalkNumsName,'chosen')
        
        % Adjust the indices in the new good stalknums variable
        load(GoodStalkNumsName,'stalknums');
        goodstalknums = stalknums;
        load(AllEllipseName,'A');
        Nstalks = length(A);
        [newgoodstalknums] = shift_indices(problem_indices,Nstalks,goodstalknums);
        
        NewGoodStalks = strcat(output_prefix,'_NewGoodStalks.mat');
        save(NewGoodStalks,'newgoodstalknums');
        
        % Run PCA
        PCA_ellipse_fits(AllGoodEllipseFits,NEPCName);

        % Create the Abaqus Python scripts
        create_cases_shifted(NEPCName,ChosenGoodEllipseFits,NewGoodStalks,problem_indices,5,material_method,group,MaterialsName);
        
    else
        % Run PCA
        PCA_ellipse_fits(AllGoodEllipseFits,NEPCName);

        % Create the Abaqus Python scripts
        create_cases(NEPCName,ChosenEllipseName,ChooseSectionsName,problem_indices,5,material_method,group,MaterialsName);
    end
    
    
end

set(0,'DefaultFigureWindowStyle','normal');

end









%% Localizing all functions used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ChooseSections(method,stalknums,dist,Table,error_indices,npoints,SaveName)
% ChooseSections.m: Determine the cross-sections to compile, which is
% determined by a method

% range: For Stalk_Table, this must be a row vector of two integer values
% from 1 to 980.

allrows = size(Table,1);

switch method
    % Choose a number of cross-sections that are all at the same distance
    % from the node
    case 'samedist'
        indices = zeros(1,length(stalknums));
%         dist = input('Choose approximate slice distance to use: ');
        
        % Step through stalks of interest (defined by range values)
        stalk = 1;
        for i = stalknums
            A = (Table.StkNum == i);
            
            % Get index of first row that is part of the current stalk
            for j = 1:length(A)
                if A(1) == 1
                    idx_first = 1;
                    break
                elseif A(j) == 1 && A(j-1) == 0
                    idx_first = j;
                end
            end
            
            % Get index of last row that is part of the current stalk
            for j = 2:length(A)
                if A(j) == 0 && A(j-1) == 1
                    idx_last = j-1;
                    break
                end
            end
            
            % Iterate through the indices that are part of the current
            % stalk and determine the index that is closest to the desired
            % value
            diffs = NaN(size(Table.StkNum));            
            for j = idx_first:idx_last
                diffs(j) = dist - Table.SlP(j);
            end
            
            % Get index of closest slice
            [~,closestIndex] = min(abs(diffs));
            
            % Add the index to indices output
            indices(stalk) = closestIndex;
            stalk = stalk + 1;
        end
        
        % Create table from indices for later reference
        selectedTable = Table(indices,:);
        
        % Get rid of any cross-sections that were chosen that are also
        % listed in error_indices
        for i = 1:length(error_indices)
            if ismember(error_indices(i),selectedTable.StkNum)
                row = find(selectedTable.StkNum == error_indices(i));
                selectedTable(row,:) = [];
            end
        end
        
        % Also get rid of any cross-sections that are made up of zeros
        % (some sort of detection error)
        deletesections = [];
        for i = 1:size(selectedTable,1)
            exterior_X = cell2mat(selectedTable.Ext_X(i));
            exterior_Y = cell2mat(selectedTable.Ext_Y(i));
            exterior_Rho = cell2mat(selectedTable.Ext_Rho(i));
            interior_X = cell2mat(selectedTable.Int_X(i));
            interior_Y = cell2mat(selectedTable.Int_Y(i));
            interior_Rho = cell2mat(selectedTable.Int_Rho(i));

            if (all(exterior_X == 0) || all(exterior_Y == 0) ||...
                    all(exterior_Rho == 0) || all(interior_X == 0) ||...
                    all(interior_Y == 0) || all(interior_Rho == 0))
                    
                deletesections = [deletesections, i];    
                msg = sprintf('Section %d deleted',i);
                disp(msg);
            end
        end
        
        %GET RID OF THE BAD SECTIONS HERE BY INDEXING
        selectedTable(deletesections,:) = [];
        
        % Save compiled slices in arrays for downstream use
        ext_X =     makearray(selectedTable,'Ext_X',npoints);
        ext_Y =     makearray(selectedTable,'Ext_Y',npoints);
        int_X =     makearray(selectedTable,'Int_X',npoints);
        int_Y =     makearray(selectedTable,'Int_Y',npoints);
        ext_T =     makearray(selectedTable,'Ext_T',npoints);
        ext_Rho =   makearray(selectedTable,'Ext_Rho',npoints);
        int_T =     makearray(selectedTable,'Int_T',npoints);
        int_Rho =   makearray(selectedTable,'Int_Rho',npoints);
        avg_rind_thick = selectedTable.rind_t;
        
        % Output all variables into mat file
        FolderName = pwd;
        SaveFile = fullfile(FolderName, SaveName);
        save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_T','ext_Rho',...
            'int_T','int_Rho','avg_rind_thick','indices','selectedTable','npoints');
        
        
        
    case 'wholestalk'
        % Choose a range of stalk numbers, and all the slices from each of
        % the chosen stalks will be chosen
        
        
        
    case 'all'
        % Chooses every slice and converts it into an array format for
        % working with more easily.
        
        % Save compiled slices in arrays for downstream use
        ext_X =     makearray(Table,'Ext_X',npoints);
        ext_Y =     makearray(Table,'Ext_Y',npoints);
        int_X =     makearray(Table,'Int_X',npoints);
        int_Y =     makearray(Table,'Int_Y',npoints);
        ext_T =     makearray(Table,'Ext_T',npoints);
        ext_Rho =   makearray(Table,'Ext_Rho',npoints);
        int_T =     makearray(Table,'Int_T',npoints);
        int_Rho =   makearray(Table,'Int_Rho',npoints);
        avg_rind_thick = Table.rind_t;
        
        % Output all variables into mat file
        FolderName = pwd;
        SaveFile = fullfile(FolderName, SaveName);
        save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_T','ext_Rho',...
            'int_T','int_Rho','avg_rind_thick','npoints');
        
    otherwise
        disp('Unknown method.');
end


end


function [exported_array] = makearray(Table,Variable,npoints)
    % Make an array of the chosen variable in the format required to work
    % with the downstream process in ellipse_fittingV1.m and
    % PCA_ellipse_fits.m
    num_slices = size(Table,1);
    exported_array = NaN(npoints,num_slices);
    
    if iscell(Table.(Variable))
        for i = 1:num_slices
            exported_array(:,i) = cell2mat(Table.(Variable)(i));
        end
    else
        for i = 1:num_slices
            exported_array(i,:) = Table.(Variable)(i);
        end
    end
        

end


function find_flip_notches(ChooseSectionsOutput,SaveName)
% FILENAME: find_flip_notches.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: Manually identify the cross-sections that need to be flipped 180
% degrees so the notch is on the left side. Cycles through the data
% contained in a .mat file output from ChooseSections.m, and allows the
% user to tag cross-sections that need to be flipped by entering "1" and
% then Enter. All other sections can be skipped through by hitting Enter.
% 
% 
% INPUTS:
%       ChooseSectionsOutput - A .mat file produced by ChooseSections.m
%       
% OUTPUTS:
%       flip_sections - A vector of 1s and empties (or non-1s) that must be
%       fed into flip_notches.m.
%
%
% NOTES: The reason find_flip_notches.m and flip_notches.m are not combined
% into a single function is so the user has the opportunity to correct for
% miskeyed cross-sections. flip_sections can be manually edited before
% being fed into flip_notches in case the user double-typed or had some
% other problem, since find_flip_notches.m only marches forward through the
% cross-sections, with no opportunity to make corrections during the
% process.
% 
% 
% VERSION HISTORY:
% V1 - Made into a function that works with the updated process flow
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load(ChooseSectionsOutput);

R_ext = ext_Rho;
R_int = int_Rho;
T = ext_T(:,1)';

N = size(R_ext,2);
flip_sections = zeros(N,1);

for i = 1:N
    % Plot each cross section to see if it needs to be flipped 180 degrees
    polarplot(T,R_ext(:,i));
    i
    s = input('Enter 1 if cross section needs to flip: ');
    if isempty(s)
        s = 0;
    end
    flip_sections(i) = s;
%     pause(); 
end
close;

% Save data as mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, SaveName);
save(SaveFile,'flip_sections');


end


function flip_notches(Flip_Indices,ChooseSectionsOutput,FlippedOutputName)
% FILENAME: find_flip_notches.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: 
% 
% 
% INPUTS:
%       flip_sections - The vector of 1s and empties that identifies the
%       cross-sections needing to be flipped 180 degrees
%
%       ChooseSectionsOutput - A .mat file produced by ChooseSections.m
%
%       SaveName - The name for the output .mat file. Make sure to end the
%       name with FLIPPED for consistency.
%       
% OUTPUTS:
%       
%
% NOTES:
%       Make sure to name the output .mat file with a FLIPPED extension
%       for consistency
% 
% 
% VERSION HISTORY:
% V1 - Made into a function that works with the updated process flow
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load(ChooseSectionsOutput);
load(Flip_Indices);

ext_rho = [];
int_rho = [];
% ext_X = ext_X;
% ext_Y = ext_Y;
% int_X = int_X;
% int_Y = int_Y;
T = ext_T(:,1)';
N = size(ext_X,2);

flipped = []; % Hold the indices of the flipped sections
for i = 1:N
    if flip_sections(i) == 1 
        flipped = [flipped; i];
        
        % Rotate and reorder external and internal points
        [~, ~, ~, ~, ~, ~, xp_ext, yp_ext, ~, ~] = reorder_V2(ext_X(:,i), ext_Y(:,i), pi);
        [~, ~, ~, ~, ~, ~, xp_int, yp_int, ~, ~] = reorder_V2(int_X(:,i), int_Y(:,i), pi);
        
        [~, ~, x_ext, y_ext, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_ext, yp_ext, 0);
        [~, ~, x_int, y_int, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_int, yp_int, 0);
        
        % Redefine the appropriate row in the main XY data
        ext_X(:,i) = x_ext;
        ext_Y(:,i) = y_ext;
        int_X(:,i) = x_int;
        int_Y(:,i) = y_int;
    end
end

% Convert data from Cartesian to polar
size(T)
for i = 1:N
    for j = 1:length(T)
        ext_Rho(j,i) = sqrt(ext_X(j,i)^2 + ext_Y(j,i)^2);
        int_Rho(j,i) = sqrt(int_X(j,i)^2 + int_Y(j,i)^2);
    end
end

% Transpose rho arrays so they are the same orientation as the other
% variables
% ext_rho = ext_rho';
% int_rho = int_rho';

flippedTable = selectedTable;

for i = 1:N
    flippedTable.Ext_X{i} = ext_X(:,i);
    flippedTable.Ext_Y{i} = ext_Y(:,i);
%     flippedTable.Ext_T{i} = ext_T(:,i);
    flippedTable.Ext_Rho{i} = ext_Rho(:,i);
    flippedTable.Int_X{i} = int_X(:,i);
    flippedTable.Int_Y{i} = int_Y(:,i);
%     flippedTable.Int_T{i} = int_T(:,i);
    flippedTable.Int_Rho{i} = int_Rho(:,i);
    
end

% Save data as mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, FlippedOutputName);
save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_Rho','int_Rho','ext_T',...
    'int_T','avg_rind_thick','flip_sections','flippedTable','npoints');

end


function ellipse_fitting_V2(FileName,SaveName)
% Load a mat file that has exterior XY data and avgrindthickness data.
% Cycle through cross sections and select the angular range that
% contains the notch so it's ignored during ellipse fitting. Then save the
% ellipses in a mat file. Take the difference between the interior and
% exterior points and their ellipse approximations and save those as well
% for later PCA.

close all

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'avg_rind_thick','ext_Rho','ext_T','ext_X','ext_Y',...
    'int_Rho','int_T','int_X','int_Y','npoints');

% Make copies of original data to work with
R_ext = ext_Rho;
R_int = int_Rho;
T = (ext_T(:,1));
X_ext = ext_X;
Y_ext = ext_Y;
X_int = int_X;
Y_int = int_Y;

N = length(avg_rind_thick);

% Instantiate output variables
A = zeros(N,1);
B = zeros(N,1);
ALPHA = zeros(N,1);

ELLIPSE_XY = zeros(N,npoints,2);
ELLIPSE_CENTERS = zeros(N,2);
ELLIPSE_T = zeros(N,npoints);
ELLIPSE_R_ext = zeros(N,npoints);
ELLIPSE_R_int = zeros(N,npoints);

DIFF_R_ext = zeros(N,npoints);
DIFF_R_int = zeros(N,npoints);
R_ext = zeros(N,npoints);
R_int = zeros(N,npoints);

AVG_RIND_T = zeros(N,1);



min_angle = 135;
min_angle = min_angle*(pi/180);     % Convert angle to radians
max_angle = 225;
max_angle = max_angle*(pi/180);

for i = 1:N
    prev_alpha = 0;
%     i
    % Define the notch range
    for j = 1:npoints
        if T(j) > min_angle
            min_index = j-1;
            break
        end
    end
    
    for j = 1:npoints
        if T(j) > max_angle
            max_index = j;
            break
        end
    end
    
    % Cut out the notch from the XY data
    window = [linspace(1,min_index,min_index),linspace(max_index,npoints,(npoints-max_index + 1))];
    xcut = X_ext(window,i);
    ycut = Y_ext(window,i);
%     if i == 892
%         xcut
%         ycut
%     end

    % Fit an ellipse to the data with the notch removed
    [alpha, major, minor, xbar_e, ybar_e, X_ellipse, Y_ellipse] = fit_ellipse_R4( xcut, ycut, npoints, prev_alpha, gca );
    
    % Save ellipse center shift
    ELLIPSE_CENTERS(i,:) = [mean(X_ellipse), mean(Y_ellipse)];
    
    % Shift XY data to be centered at the origin before converting to polar
    X_ellipse_shift = X_ellipse - mean(X_ellipse);
    Y_ellipse_shift = Y_ellipse - mean(Y_ellipse);
    X_ext_shift = X_ext(:,i) - mean(X_ellipse);
    Y_ext_shift = Y_ext(:,i) - mean(Y_ellipse);
    X_int_shift = X_int(:,i) - mean(X_ellipse);
    Y_int_shift = Y_int(:,i) - mean(Y_ellipse);
    
    % Reorder indices to start at 0 degrees
    [X_ellipse_shift, Y_ellipse_shift, ~, ~, ~, ~, ~, ~, ~, ~] = reorder_V2(X_ellipse_shift, Y_ellipse_shift, 0);
    [X_ext_shift, Y_ext_shift, ~, ~, ~, ~, ~, ~, ~, ~] = reorder_V2(X_ext_shift, Y_ext_shift, 0);
    [X_int_shift, Y_int_shift, ~, ~, ~, ~, ~, ~, ~, ~] = reorder_V2(X_int_shift, Y_int_shift, 0);    
    
    % Convert X_ellipse and Y_ellipse to polar coordinates
    theta = 0:2*pi/npoints:2*pi;
    theta = theta(1:end-1);
    [thetatemp_ellipse, ext_rho_ellipse] = cart2pol(X_ellipse_shift,Y_ellipse_shift);
    thetatemp_ellipse = wrapTo2Pi(thetatemp_ellipse);   % Make all the negative pi values sit on a 0-2*pi system
    
    % Remove duplicate values for interpolation
    [C,ia,~] = unique(thetatemp_ellipse,'stable');
    thetatemp_ellipse = C;
    ext_rho_ellipse = ext_rho_ellipse(ia);
    
    % Convert X_ext_shift and Y_ext_shift to polar for resampling
    [thetatemp_ext, ext_rho] = cart2pol(X_ext_shift,Y_ext_shift);
    thetatemp_ext = wrapTo2Pi(thetatemp_ext);   % Make all the negative pi values sit on a 0-2*pi system
    
    % Remove duplicate values for interpolation
    [C,ia,~] = unique(thetatemp_ext,'stable');
    thetatemp_ext = C;
    ext_rho = ext_rho(ia);
    
    % Convert X_int_shift and Y_int_shift to polar for resampling
    [thetatemp_int, int_rho] = cart2pol(X_int_shift,Y_int_shift);
    thetatemp_int = wrapTo2Pi(thetatemp_int);   % Make all the negative pi values sit on a 0-2*pi system
    
    % Remove duplicate values for interpolation
    [C,ia,~] = unique(thetatemp_int,'stable');
    thetatemp_int = C;
    int_rho = int_rho(ia);
    
    
    
    
    % Interpolate to get new rho and theta points that are regularly spaced
    ext_rho_ellipse_interp = interp1(thetatemp_ellipse,ext_rho_ellipse,theta,'pchip','extrap'); 
    ext_rho_ellipse = ext_rho_ellipse_interp;
    ext_rho_interp = interp1(thetatemp_ext,ext_rho,theta,'pchip','extrap'); 
    ext_rho = ext_rho_interp;
    int_rho_interp = interp1(thetatemp_int,int_rho,theta,'pchip','extrap'); 
    int_rho = int_rho_interp;
    
    % Get interior ellipse fit points, based on constant rind thickness
    % assumption
    int_rho_ellipse = ext_rho_ellipse - avg_rind_thick(i);
    
%     % Plot in polar coordinates to check results
%     polarplot(theta,ext_rho,'.','LineWidth',2);
%     hold on
%     polarplot(theta,int_rho,'.','LineWidth',2);
% %     polarplot(theta,ext_rho_ellipse,'.','LineWidth',2)
%     pause();
%     close;
    
    
    
    A(i) = major;
    B(i) = minor;
    ALPHA(i) = alpha;

    ELLIPSE_XY(i,:,1) = X_ellipse;
    ELLIPSE_XY(i,:,2) = Y_ellipse;
    ELLIPSE_CENTERS(i,1) = mean(X_ellipse);
    ELLIPSE_CENTERS(i,2) = mean(Y_ellipse);
    ELLIPSE_T(i,:) = theta;
    ELLIPSE_R_ext(i,:) = ext_rho_ellipse;
%     ELLIPSE_R_int(i,:) = ext_rho_ellipse - avg_rind_thick(i);
    ELLIPSE_R_int(i,:) = rpts(npoints,ELLIPSE_T(i,:),(A(i) - 2*avg_rind_thick(i)),(B(i) - 2*avg_rind_thick(i)));
    R_ext(i,:) = ext_rho;
    R_int(i,:) = int_rho;
    
    % Get difference between the ellipse and the real data (if the ellipse
    % overestimates, then the value of DIFF will be positive)
    DIFF_R_ext(i,:) = ext_rho_ellipse - ext_rho;
    DIFF_R_int(i,:) = int_rho_ellipse - int_rho;
    
    AVG_RIND_T(i) = avg_rind_thick(i);
    
    RIND_ELLIPSE_DIFF(i,:) = ELLIPSE_R_ext(i,:) - R_int(i,:);
    
end

% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int',...
    'ELLIPSE_CENTERS','DIFF_R_ext','DIFF_R_int','R_ext','R_int','AVG_RIND_T');

end


function [r] = rpts(N,theta,dmaj,dmin)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end


function PCA_ellipse_fits(FileName,SaveName)
% USE THIS FUNCTION ON Ellipse_fits_bottom1.mat or Ellipse_fits_top1.mat
% (uses the difference between the ellipse and the real data)

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'DIFF_R_ext','DIFF_R_int','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int');

% Perform PCA. 'Centered' option must be false to allow for reverse
% engineering of the original data
[ext_rhoPCAs, ext_rhocoeffs, ext_rhoPCA_variances, ext_rhotstat, ext_rhoexplained, ext_rhovarMeans] = pca(DIFF_R_ext,'Centered',false);
% [int_rhoPCAs, int_rhocoeffs, int_rhoPCA_variances, int_rhotstat, int_rhoexplained, int_rhovarMeans] = pca(DIFF_R_int,'Centered',false);

ext_rhoexplained_tot = zeros(size(ext_rhoexplained));
% int_rhoexplained_tot = zeros(size(ext_rhoexplained));
for i = 1:length(ext_rhoexplained_tot)
    ext_rhoexplained_tot(i) = sum(ext_rhoexplained(1:i));
%     int_rhoexplained_tot(i) = sum(int_rhoexplained(1:i));
end

figure(1);
plot(ext_rhoexplained_tot(:,1),'-*');
title('Exterior Non-Ellipse PCs');
xlabel('# of PCs');
ylabel('% Variance Explained');

% figure(2);
% plot(int_rhoexplained_tot(:,1),'-*');
% title('Interior Non-Ellipse PCs');
% xlabel('# of PCs');
% ylabel('% Variance Explained');

ELLIPSE_T = ELLIPSE_T';
theta = ELLIPSE_T(:,1);
ELLIPSE_T = ELLIPSE_T';

figure(3);
polarplot(theta,ext_rhoPCAs(:,1));
hold on
polarplot(theta,ext_rhoPCAs(:,2));
polarplot(theta,ext_rhoPCAs(:,3));
polarplot(theta,ext_rhoPCAs(:,4));
polarplot(theta,ext_rhoPCAs(:,5));
title('Exterior Rho Principal Components');
legend('PC1','PC2','PC3','PC4','PC5');

% figure(4);
% polarplot(theta,int_rhoPCAs(:,1));
% hold on
% polarplot(theta,int_rhoPCAs(:,2));
% polarplot(theta,int_rhoPCAs(:,3));
% polarplot(theta,int_rhoPCAs(:,4));
% polarplot(theta,int_rhoPCAs(:,5));
% title('Interior Rho Principal Components');
% legend('PC1','PC2','PC3','PC4','PC5');




% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'ELLIPSE_T','ELLIPSE_R_ext','ext_rhocoeffs',...
    'ext_rhoPCAs','ext_rhoexplained','ext_rhovarMeans');


end


function remove_problem_ellipses(OriginalEllipseFits,stalknums,problem_indices,GoodEllipseFits,GoodStalks,All)
load(OriginalEllipseFits);
% load(Stalks,'stalknums');

% Remove the problem ellipses from the ellipse fit data
A(problem_indices) = [];
AVG_RIND_T(problem_indices) = [];
B(problem_indices) = [];
DIFF_R_ext(problem_indices,:) = [];
DIFF_R_int(problem_indices,:) = [];
ELLIPSE_CENTERS(problem_indices,:) = [];
ELLIPSE_R_ext(problem_indices,:) = [];
ELLIPSE_R_int(problem_indices,:) = [];
ELLIPSE_T(problem_indices,:) = [];
ELLIPSE_XY(problem_indices,:,:) = [];
R_ext(problem_indices,:) = [];
R_int(problem_indices,:) = [];

% If removing translated_problem_indices from the chosen set of stalks
% only:
if ~strcmp(All,'all')
    if ~isempty(problem_indices)
        stalknums(problem_indices) = []; % This line should only be executed if problem_indices is translated_problem_indices
    end
end

% Save the final data in a new mat file
FolderName = pwd;
SaveFile       = fullfile(FolderName, GoodEllipseFits);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int',...
    'ELLIPSE_CENTERS','DIFF_R_ext','DIFF_R_int','R_ext','R_int','AVG_RIND_T','problem_indices');
SaveFile = fullfile(FolderName,GoodStalks);
save(SaveFile,'stalknums');

end

function create_cases_shifted(NEPCdata,GoodEllipseData,NewGoodStalks,problem_indices,numNEPCs,material_method,group,SaveName)
    % create_cases.m: Calculate the necessary information to include in the
    % Python scripts
    
    load(NEPCdata);
    load(GoodEllipseData);
    load(NewGoodStalks,'newgoodstalknums'); % Load in the shifted stalknums

    N = size(ELLIPSE_T,1);
    MaterialProps = zeros(N,(2 + 2*numNEPCs - 1),2);
    
%     if strcmp(cpu,'multi')
%         write_Python_template_mult_processor;
%     elseif license == 1
%         write_Python_template;  % Create Template cell array that can be copied and used to make individualized Python scripts
%     elseif license == 2
%         write_Python_template2;
%     else
%         error('Choose a valid cpu/license combination');
%     end

    write_Python_template3;

    %% Create all geometry cases for a given cross section
    % Step through the cross sections
%     stalknums = selectedTable.StkNum; % NEED A TRANSLATED stalks VARIABLE SO IT DOESN'T CHOOSE CROSS-SECTIONS OUTSIDE
    
%     % Remove the stalk numbers that had ellipse fit problems REDUNDANT
%     if ~isempty(problem_indices)
%         stalknums(stalknums==problem_indices) = [];
%     end
    
    for i = 1:N
        GROUP = sprintf('%d',group); % Group number
        ID = sprintf('%d',newgoodstalknums(i)); % Cross-section number

        %% Real cross section (case 0)
        case_num = 0; % increment this for each case within each cross section
        Script = Template;
        [Erind,Epith] = get_materials(material_method);
        make_case(case_num,i,ID,GROUP,R_ext,R_int,ELLIPSE_T,Script,Erind,Epith);
        MaterialProps(i,case_num+1,1) = Erind;
        MaterialProps(i,case_num+1,2) = Epith;

        %% Pure ellipse fit (case 1)
        case_num = case_num + 1;
        Script = Template; % Reset the script template
        [Erind,Epith] = get_materials(material_method);
        make_case(case_num,i,ID,GROUP,ELLIPSE_R_ext,ELLIPSE_R_int,ELLIPSE_T,Script,Erind,Epith);
        MaterialProps(i,case_num+1,1) = Erind;
        MaterialProps(i,case_num+1,2) = Epith;

        %% Combined NEPC cases
        for j = 1:numNEPCs
            case_num = case_num + 1;
            Script = Template; % Reset the script template

            % Calculate the cases with NEPCs cumulatively added into the
            % ellipse fit
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
%             NEPC_int = zeros(1,size(int_rhoPCAs,1));
            for k = 1:j
                % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
                NEPC_ext = NEPC_ext + ext_rhocoeffs(newgoodstalknums(i),k)*ext_rhoPCAs(:,k)';
%                 NEPC_int = NEPC_int + int_rhocoeffs(i,k)*int_rhoPCAs(:,k)';
            end

            Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
            Rnew_int = Rnew_ext - AVG_RIND_T(i);
            
            [Erind,Epith] = get_materials(material_method);
            make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,ELLIPSE_T,Script,Erind,Epith);
            MaterialProps(i,case_num+1,1) = Erind;
            MaterialProps(i,case_num+1,2) = Epith;

        end


        %% Remaining individual NEPC cases
        for j = 2:numNEPCs
            case_num = case_num + 1;
            Script = Template; % Reset the script template

            % Add the current NEPC to the ellipse in polar coordinates
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
%             NEPC_int = zeros(1,size(int_rhoPCAs,1));
            NEPC_ext = ext_rhocoeffs(newgoodstalknums(i),j)*ext_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
%             NEPC_int = int_rhocoeffs(i,j)*int_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
            Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
%             Rnew_int = Rnew_ext - AVG_RIND_T(i);

            [Erind,Epith] = get_materials(material_method);
            make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,ELLIPSE_T,Script,Erind,Epith);
            MaterialProps(i,case_num+1,1) = Erind;
            MaterialProps(i,case_num+1,2) = Epith;

        end

    end
    
    % Save the final data in a new mat file
    FolderName = pwd;
    SaveFile = fullfile(FolderName, SaveName);
    save(SaveFile,'MaterialProps');
    
    

end


function create_cases(NEPCdata,GoodEllipseData,ChosenSectionsData,problem_indices,numNEPCs,material_method,group,SaveName)
    % create_cases.m: Calculate the necessary information to include in the
    % Python scripts
    
    load(NEPCdata);
    load(GoodEllipseData);
    load(ChosenSectionsData);
%     load(NewGoodStalks,'newgoodstalknums'); % Load in the shifted stalknums

    N = size(ELLIPSE_T,1);
    MaterialProps = zeros(N,(2 + 2*numNEPCs - 1),2);
    
%     if strcmp(cpu,'multi')
%         write_Python_template_mult_processor;
%     elseif license == 1
%         write_Python_template;  % Create Template cell array that can be copied and used to make individualized Python scripts
%     elseif license == 2
%         write_Python_template2;
%     else
%         error('Choose a valid cpu/license combination');
%     end

    write_Python_template3;

    %% Create all geometry cases for a given cross section
    % Step through the cross sections
%     stalknums = selectedTable.StkNum; % NEED A TRANSLATED stalks VARIABLE SO IT DOESN'T CHOOSE CROSS-SECTIONS OUTSIDE
    
    % Remove the stalk numbers that had ellipse fit problems REDUNDANT
    if ~isempty(problem_indices)
        stalknums(stalknums==problem_indices) = [];
    end
    
    for i = 1:N
        GROUP = sprintf('%d',group); % Group number
        ID = sprintf('%d',stalknums(i)); % Cross-section number

        %% Real cross section (case 0)
        case_num = 0; % increment this for each case within each cross section
        Script = Template;
        [Erind,Epith] = get_materials(material_method);
        make_case(case_num,i,ID,GROUP,R_ext,R_int,ELLIPSE_T,Script,Erind,Epith);
        MaterialProps(i,case_num+1,1) = Erind;
        MaterialProps(i,case_num+1,2) = Epith;

        %% Pure ellipse fit (case 1)
        case_num = case_num + 1;
        Script = Template; % Reset the script template
        [Erind,Epith] = get_materials(material_method);
        make_case(case_num,i,ID,GROUP,ELLIPSE_R_ext,ELLIPSE_R_int,ELLIPSE_T,Script,Erind,Epith);
        MaterialProps(i,case_num+1,1) = Erind;
        MaterialProps(i,case_num+1,2) = Epith;

        %% Combined NEPC cases
        for j = 1:numNEPCs
            case_num = case_num + 1;
            Script = Template; % Reset the script template

            % Calculate the cases with NEPCs cumulatively added into the
            % ellipse fit
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
%             NEPC_int = zeros(1,size(int_rhoPCAs,1));
            for k = 1:j
                % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
                NEPC_ext = NEPC_ext + ext_rhocoeffs(newgoodstalknums(i),k)*ext_rhoPCAs(:,k)';
%                 NEPC_int = NEPC_int + int_rhocoeffs(i,k)*int_rhoPCAs(:,k)';
            end

            Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
            Rnew_int = Rnew_ext - AVG_RIND_T(i);
            
            [Erind,Epith] = get_materials(material_method);
            make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,ELLIPSE_T,Script,Erind,Epith);
            MaterialProps(i,case_num+1,1) = Erind;
            MaterialProps(i,case_num+1,2) = Epith;

        end


        %% Remaining individual NEPC cases
        for j = 2:numNEPCs
            case_num = case_num + 1;
            Script = Template; % Reset the script template

            % Add the current NEPC to the ellipse in polar coordinates
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
%             NEPC_int = zeros(1,size(int_rhoPCAs,1));
            NEPC_ext = ext_rhocoeffs(stalknums(i),j)*ext_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
%             NEPC_int = int_rhocoeffs(i,j)*int_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
            Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
%             Rnew_int = Rnew_ext - AVG_RIND_T(i);

            [Erind,Epith] = get_materials(material_method);
            make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,ELLIPSE_T,Script,Erind,Epith);
            MaterialProps(i,case_num+1,1) = Erind;
            MaterialProps(i,case_num+1,2) = Epith;

        end

    end
    
    % Save the final data in a new mat file
    FolderName = pwd;
    SaveFile = fullfile(FolderName, SaveName);
    save(SaveFile,'MaterialProps');
    
    

end

function make_case(case_num,i,ID,GROUP,R_ext,R_int,T,Script,Erind,Epith)
    CASE = sprintf('%d',case_num);
    jobname = strcat('''Group_',GROUP,'_','Section_',ID,'_',CASE,'''');
    scriptname = strcat('Group_',GROUP,'_','Section_',ID,'_',CASE,'.py');
    
    % Convert data to Cartesian coordinates (read in as row vectors)
    if size(R_ext,1) > 1
        X_ext = R_ext(i,:).*cos(T(i,:));
        Y_ext = R_ext(i,:).*sin(T(i,:));
        X_int = R_int(i,:).*cos(T(i,:));
        Y_int = R_int(i,:).*sin(T(i,:));
    else
        X_ext = R_ext(1,:).*cos(T(1,:));
        Y_ext = R_ext(1,:).*sin(T(1,:));
        X_int = R_int(1,:).*cos(T(1,:));
        Y_int = R_int(1,:).*sin(T(1,:));
    end

     % Scale units to micrometers from millimeters
    X_ext = 1000*X_ext;
    Y_ext = 1000*Y_ext;
    X_int = 1000*X_int;
    Y_int = 1000*Y_int;
    
    
    % Transpose data and combine xy
    section_ext = [X_ext', Y_ext'];
    section_int = [X_int', Y_int'];

    % Repeat the last points to close the loop
    section_ext = [section_ext; section_ext(1,:)];
    section_int = [section_int; section_int(1,:)];

    % Get the reference point values in Cartesian coordinates for
    % reference points closest to 90 and 270 degrees
    diffs90 = NaN(1,size(T,2));
    diffs270 = NaN(1,size(T,2));
    for j = 1:length(T(1,:))
        diffs90(j) = pi/2 - T(1,j);
        diffs270(j) = 3*pi/2 - T(1,j);
    end
    
    [~,ind90] = min(abs(diffs90));
    [~,ind270] = min(abs(diffs270));
    
    RP1X = sprintf('%0.5g',X_ext(ind90));
    RP1Y = sprintf('%0.5g',Y_ext(ind90));
    RP2X = sprintf('%0.5g',X_ext(ind270));
    RP2Y = sprintf('%0.5g',Y_ext(ind270));

    % Write the spline points and save as a string
    S = size(section_ext);
    len = S(1);
    outer_spline = writespline_V2(len,section_ext);
    inner_spline = writespline_V2(len,section_int);
    
%     % Calculate the random material properties from a normal distribution.
%     % Bound with 95% confidence interval, calculated from transverse
%     % material properties used in another paper.
%     Erind_mean = 8.0747e-04;
%     Erind_stdev = 3.3517e-04;
%     Erind_95 = [6.7414e-04 9.4081e-04];
%     Epith_mean = 2.5976e-05;
%     Epith_stdev = 1.0303e-05;
%     Epith_95 = [2.1878e-05 3.0075e-05];
%     
%     % Generate Erind from normal distribution
%     while 1
%         Erind = normrnd(Erind_mean,Erind_stdev);
%         if Erind >= Erind_95(1) && Erind <= Erind_95(2)
%             break
%         end
%     end
%     
%     % Generate Epith from normal distribution
%     while 1
%         Epith = normrnd(Epith_mean,Epith_stdev);
%         if Epith >= Epith_95(1) && Epith <= Epith_95(2)
%             break
%         end
%     end
    
    rindE = sprintf('%0.5g',Erind);
    pithE = sprintf('%0.5g',Epith);

    % Insert the case-specific values into the appropriate parts of the
    % Python script template (must be strings)
    Script(17,1) = strcat(Script(17,1),jobname);
    Script(19,1) = strcat(Script(19,1),GROUP);
    Script(23,1) = strcat(Script(23,1),ID);
    Script(25,1) = strcat(Script(25,1),CASE);
    Script(33,1) = strcat(Script(33,1),rindE);
    Script(35,1) = strcat(Script(35,1),pithE);
    Script(37,1) = strcat(Script(37,1),RP1X);
    Script(39,1) = strcat(Script(39,1),RP1Y);
    Script(41,1) = strcat(Script(41,1),RP2X);
    Script(43,1) = strcat(Script(43,1),RP2Y);
    Script(63,1) = strcat(Script(63,1),outer_spline);
    Script(86,1) = strcat(Script(86,1),inner_spline);
    
    % Write Python script from the cell array
    filePh = fopen(scriptname,'w');
    fprintf(filePh,'%s\n',Script{:});
    fclose(filePh);
    
end

function [spline] = writespline_V2(len,data)
    %define empty spline and number of x-y points
    spline = '';

    %run through 1-column arrays of the x and y data points for the spline, and add to the end of the string with the correct formatting
    for i = 1:len 
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), '); 
    end
end

function [xy_columns] = convert_to_xy(R,theta)
    N = length(theta);
    xy_columns = zeros(N,2);
    for i = 1:N
        xy_columns(i,1) = R(i)*cos(theta(i));
        xy_columns(i,2) = R(i)*sin(theta(i));
    end
end

function [Erind,Epith] = get_materials(method)
% Calculate the random material properties from a normal distribution.
    % Bound with 95% confidence interval, calculated from transverse
    % material properties used in another paper.
    Erind_mean = 8.0747e-04;
    Erind_stdev = 3.3517e-04;
    Erind_95 = [6.7414e-04 9.4081e-04];
    Epith_mean = 2.5976e-05;
    Epith_stdev = 1.0303e-05;
    Epith_95 = [2.1878e-05 3.0075e-05];
    ratio_mean = 0.0372;
    ratio_stdev = 0.0180;
    ratio_95 = [0.0300 0.0444];
    
    switch method
        case 'random'
            % Generate Erind from normal distribution
            while 1
                Erind = normrnd(Erind_mean,Erind_stdev);
                if Erind >= Erind_95(1) && Erind <= Erind_95(2)
                    break
                end
            end

            % Generate Epith from normal distribution
            while 1
                Epith = normrnd(Epith_mean,Epith_stdev);
                if Epith >= Epith_95(1) && Epith <= Epith_95(2)
                    break
                end 
            end
            
            
    %     % Generate Epith from normal distribution of pith/rind ratios
    %     while 1
    %         ratio = normrnd(ratio_mean,ratio_stdev);
    %         if ratio >= ratio_95(1) && ratio <= ratio_95(2)
    %             break
    %         end
    %     end
    %     Epith = ratio*Erind;

    
        case 'min'
            Erind = Erind_95(1);
            Epith = Epith_95(1);
            
        case 'max'
            Erind = Erind_95(2);
            Epith = Epith_95(2);
            
        case 'minpith'
            Erind = Erind_mean;
            Epith = Epith_95(1);
            
        case 'maxpith'
            Erind = Erind_mean;
            Epith = Epith_95(2);
            
        case 'minrind'
            Erind = Erind_95(1);
            Epith = Epith_mean;
            
        case 'maxrind'
            Erind = Erind_95(2);
            Epith = Epith_mean;
    
        case 'avg'
            Erind = Erind_mean;
            Epith = Epith_mean;
    end
    
end