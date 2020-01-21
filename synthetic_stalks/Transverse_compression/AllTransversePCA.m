function AllTransversePCA(slice_dists)
% FILENAME: AllTransversePCA.m
% AUTHOR: Ryan Larson
% DATE: 1/17/2020
%
% PURPOSE: 
% 
% 
% INPUTS:
%       - slice_dists: A row vector of all the slice distances to take into
%       acccount when running PCA (MAYBE MAKE IT POSSIBLE TO LATER CHOOSE
%       WHICH SLICE DISTANCES ACTUALLY GET USED IN THE PCA GOING INTO
%       TRANSVERSE_WRAPPER_V4.M)
%       
% OUTPUTS:
%       
%
%
% NOTES: 
%       
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
load StalksDCR_360pts.mat
hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

ALL_PROBLEM_INDICES = [];
ALL_DIFF_R_ext      = [];
ALL_DIFF_R_int      = [];
ALL_ELLIPSE_T       = [];
ALL_ELLIPSE_R_ext   = [];
ALL_ELLIPSE_R_int   = [];

slice_startstop = zeros(length(slice_dists),3);
slice_startstop(:,1) = slice_dists';

for slice = slice_dists
    slice
    % For each slice distance, iterate through the 980 stalks and get the
    % data to feed into the large PCA array
    
    % Set up naming for .mat files
    dist_int = num2str(floor(abs(slice)));
    deci = abs(slice) - floor(abs(slice));
    dist_deci = num2str(deci);
    dist_deci = erase(dist_deci,'0.');

    if slice > 0
        if strcmp(dist_deci,'0')
            slicepos = strcat('Above_',dist_int);
        else
            slicepos = strcat('Above_',dist_int,'_',dist_deci);
        end
    elseif slice < 0
        if strcmp(dist_deci,'0')
            slicepos = strcat('Below_',dist_int);
        else
            slicepos = strcat('Below_',dist_int,'_',dist_deci);
        end
    else
        slicepos = strcat('At_Node');
    end

    output_prefix = strcat(slicepos);
    
    % Gather all data for stalks at this slice distance
    AllSectionsName = strcat(output_prefix,'_All980.mat');
    ChooseSections('samedist',linspace(1,980,980),slice,Stalk_TableDCR,error_indices,npoints,AllSectionsName)
    load(AllSectionsName);

    
    
    %% Flip cross-sections that need adjustment (only if there isn't flip data)
    % Check to see if there is already a flip vector for the chosen distance
    % Get file name to look for
    searchslice = sprintf('%d',abs(slice));
    if slice > 0
        FlippedOutputName = strcat('Above_',searchslice,'_FLIPPED.mat');
        fstruct = dir(FlippedOutputName);
    elseif slice < 0
        FlippedOutputName = strcat('Below_',searchslice,'_FLIPPED.mat');
        fstruct = dir(FlippedOutputName);
    else
        FlippedOutputName = strcat('At_Node','_FLIPPED.mat');
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
    
    

    %% Get ellipse fits and difference data 
    %% Check to see if there is already an ellipse fit for the chosen distance
    % Get file name to look for
    searchslice = sprintf('%d',slice);
    if slice > 0
        searchname = strcat('*Above_',searchslice,'_AllGoodEllipses.mat');
        fstruct = dir(searchname);
        if isempty(fstruct)
            searchname = strcat('*Above_',searchslice,'_AllEllipses.mat');
        end

    elseif slice < 0
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
%     chosen_problem_indices = [];

    if ~isfile(searchname)
        disp('No ellipse data for this slice distance exists in the current folder. Create data now.');

        % Calculate the ellipse fits for all cross-sections at chosen distance
        AllEllipseName = strcat(output_prefix,'_AllEllipses.mat');
        ellipse_fitting_V2(FlippedOutputName,AllEllipseName);

%         % Calculate the ellipse fits for the chosen set of cross-sections
%         ChosenEllipseName = strcat(output_prefix,'_ChosenEllipses.mat');
%         ellipse_fitting_V2(ChooseSectionsName,ChosenEllipseName);

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



%         % Determine if any of the problem cases are in the selected stalks
%         problocs = ismember(problem_indices,stalknums);
%         chosen_problem_indices = problem_indices(problocs);

        % Save problem_indices for later use, after correcting the indices
        ProblemEllipses = strcat(output_prefix,'_ProblemEllipses.mat');
        save(ProblemEllipses,'problem_indices');

        problem_indices
%         chosen_problem_indices


    else
        disp('Ellipse data for the chosen slice distance has been found.');
        AllEllipseName = strcat(output_prefix,'_AllEllipses.mat');

        ProblemEllipses = strcat(output_prefix,'_ProblemEllipses.mat');
        load(ProblemEllipses,'problem_indices');

%         % Calculate the ellipse fits for the chosen set of cross-sections
%         ChosenEllipseName = strcat(output_prefix,'_ChosenEllipses.mat');
%         ellipse_fitting_V2(ChooseSectionsName,ChosenEllipseName);
    end
    
    
    %% Prepare data for PCA
    % Add the difference data and corresponding ellipse fits to separate
    % large arrays for feeding into PCA. Take note of the starting and
    % ending indices corresponding to each slice distance so data can be
    % correctly reconstructed later. Also save the error_indices that
    % correspond.
    
    load(AllEllipseName,'DIFF_R_ext','DIFF_R_int','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int');
    load(ProblemEllipses,'problem_indices');
    
    n_ALL_PROBLEM_INDICES = length(ALL_PROBLEM_INDICES);
    n_ALL_DIFF_R_ext      = size(ALL_DIFF_R_ext,1);
    n_ALL_DIFF_R_int      = size(ALL_DIFF_R_int,1);
    n_ALL_ELLIPSE_T       = size(ALL_ELLIPSE_T,1);
    n_ALL_ELLIPSE_R_ext   = size(ALL_ELLIPSE_R_ext);
    n_ALL_ELLIPSE_R_int   = size(ALL_ELLIPSE_R_int);
    
    ALL_PROBLEM_INDICES = [ALL_PROBLEM_INDICES; problem_indices];
    ALL_DIFF_R_ext      = [ALL_DIFF_R_ext; DIFF_R_ext];
    ALL_DIFF_R_int      = [ALL_DIFF_R_int; DIFF_R_int];
    ALL_ELLIPSE_T       = [ALL_ELLIPSE_T; ELLIPSE_T];
    ALL_ELLIPSE_R_ext   = [ALL_ELLIPSE_R_ext; ELLIPSE_R_ext];
    ALL_ELLIPSE_R_int   = [ALL_ELLIPSE_R_int; ELLIPSE_R_int];
    
    slice_startstop(
    
    
    % Remove error cases and adjust indices, starting at the bottom of the
    % array. Save adjusted indices.
    
    
    % Run PCA on the resulting large data set. Save this for access by
    % transverse_wrapper_V4.m
    
    
    % MAKE SURE TO ADDRESS AND CHECK THE INTERPOLATION ISSUE THAT MIGHT
    % HAVE BEEN CAUSING ALL THE ERROR_INDICES AND THE NEED TO SHIFT
    % EVERYTHING
    
end


end




%% Localized functions
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
%     int_rho_ellipse = ext_rho_ellipse - avg_rind_thick(i);
    
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
%     ELLIPSE_R_int(i,:) = ext_rho_ellipse - avg_rind_thick(i); % This way
%     gives varying rind thickness because the theta line isn't always
%     normal to the ellipse
%     ELLIPSE_R_int(i,:) = rpts(npoints,ELLIPSE_T(i,:),(A(i) - 2*avg_rind_thick(i)),(B(i) - 2*avg_rind_thick(i)));    % This way overestimates and underestimates periodically
    ELLIPSE_R_int(i,:) = normintV2(ELLIPSE_R_ext(i,:),ELLIPSE_T(i,:),avg_rind_thick(i));
    R_ext(i,:) = ext_rho;
    R_int(i,:) = int_rho;
    
    % Get difference between the ellipse and the real data (if the ellipse
    % overestimates, then the value of DIFF will be positive)
%     DIFF_R_ext(i,:) = ext_rho_ellipse - ext_rho;
%     DIFF_R_int(i,:) = int_rho_ellipse - int_rho;
    DIFF_R_ext(i,:) = ELLIPSE_R_ext(i,:) - R_ext(i,:);
    DIFF_R_int(i,:) = ELLIPSE_R_int(i,:) - R_int(i,:);
    
    AVG_RIND_T(i) = avg_rind_thick(i);
    
%     RIND_ELLIPSE_DIFF(i,:) = ELLIPSE_R_ext(i,:) - R_int(i,:);
    
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
