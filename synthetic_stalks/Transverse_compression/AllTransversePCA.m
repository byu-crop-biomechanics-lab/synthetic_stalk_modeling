function AllTransversePCA(slice_dists,SaveName)
% FILENAME: AllTransversePCA.m
% AUTHOR: Ryan Larson
% DATE: 1/17/2020
%
% PURPOSE: 
% 
% 
% INPUTS:
%       slice_dists: A row vector of all the slice distances to take into
%       acccount when running PCA (MAYBE MAKE IT POSSIBLE TO LATER CHOOSE
%       WHICH SLICE DISTANCES ACTUALLY GET USED IN THE PCA GOING INTO
%       TRANSVERSE_WRAPPER_V4.M)
%       
%       SaveName: String name to save the data to as a .mat file. Must
%       include .mat in the name.
%       
% OUTPUTS:
%       N/A
%
%
% NOTES: 
%       Other outputs:
%           - Named .mat file that contains all the ellipse fit and PCA
%           data, which will be used when running transverse_wrapper_V4.m.
%       
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
load StalksDCR_360pts_V2.mat
hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

ALL_PROBLEM_INDICES = {};
ALL_R_ext           = [];
ALL_R_int           = [];
ALL_A               = [];
ALL_B               = [];
ALL_DIFF_R_ext      = [];
ALL_DIFF_R_int      = [];
ALL_ELLIPSE_T       = [];
ALL_AVG_RIND_T      = [];
ALL_ELLIPSE_R_ext   = [];
ALL_ELLIPSE_R_int   = [];

adj_indices         = {};

slice_startstop = zeros(length(slice_dists),3);
slice_startstop(:,1) = slice_dists';

for slice = slice_dists
    slice
    close all;
    
    % For each slice distance, iterate through the 980 stalks and get the
    % data to feed into the large PCA array
    
    % Set up naming for .mat files
    dist_int = num2str(floor(abs(slice)));
    deci = abs(slice) - floor(abs(slice));
    dist_deci = num2str(deci);
    dist_deci = erase(dist_deci,'0.');
    
    % Turn the slice position into a readable string
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
    
    % Create the naming prefix for the current slice location
    output_prefix = strcat(slicepos);
    
    % Gather all data for stalks at this slice distance
    AllSectionsName = strcat(output_prefix,'_All980.mat');
    ChooseSections('samedist',linspace(1,980,980),slice,Stalk_TableDCR,error_indices,npoints,AllSectionsName)
    load(AllSectionsName);

    
    
    %% Flip cross-sections that need adjustment (only if there isn't flip data)
    % Explanation: Sometimes the downsample, center, and rotate processing
    % doesn't work exactly as intended. For principal component analysis to
    % work, all the notches need to be on the same side (the left side, in
    % the case of this code). It's easiest to do these corrections
    % manually, so this section of code allows quick visual tagging of
    % cases where the notch didn't end up on the correct side. For tagged
    % slices, it flips the cross-section by 180 degrees. There might also
    % be cases where the ellipse is vertically-oriented or where the
    % ellipse didn't do a good job of fitting. Those cases are addressed in
    % a later section of the code.    
    
    
    % Check to see if there is already flipping data for the current slice
    % location
    
    % Get a string of the slice location
    searchslice = sprintf('%d',abs(slice));
    
    % Get the appropriate filename to search for, depending on the slice
    % location
    if slice > 0
        FlippedOutputName = strcat('*Above_',searchslice,'_FLIPPED.mat');
        fstruct = dir(FlippedOutputName);
    elseif slice < 0
        FlippedOutputName = strcat('*Below_',searchslice,'_FLIPPED.mat');
        fstruct = dir(FlippedOutputName);
    else
        FlippedOutputName = strcat('*At_Node','_FLIPPED.mat');
        fstruct = dir(FlippedOutputName);        
    end
    
    % If there isn't existing flip data for the chosen slice location, then
    % create the file name.
    if isempty(fstruct)
        FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
    else
        FlippedOutputName = fstruct(1).name;
    end
    
    % FlipName will contain data about which cross-sections need to be
    % flipped, while FlippedOutputName will contain the corrected, or
    % flipped, data that will be used to perform later steps.
    FlipName = strcat(output_prefix,'_flip_sections.mat');
    
    % If there isn't existing flipped data, then create it.
    if ~isfile(FlippedOutputName)
        disp('No flip index vector exists in the current folder. Create one now.');
        
        % Manually find the cross-sections that need to be flipped 180
        % degrees (more details are in the local function section)
        find_flip_notches(AllSectionsName,FlipName)

        while 1
            % Sometimes you'll incorrectly tag a cross-section for
            % flipping, when it was already oriented correctly. This code
            % gives you a chance to make those corrections by directly
            % editing the "flip vector," or the indicators for each
            % cross-section. If you notice that you incorrectly tagged a
            % certain cross-section, note the index, which will be
            % displayed in the command window. Then the following code will
            % open the flip_sections variable, where you can change the
            % value for the incorrectly-tagged cross-section. If it should
            % be flipped, change the value to 1. If it was incorrectly
            % tagged, change the value back to 0. Follow the prompts in the
            % command window to return to the execution of the main code.
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


        % Flip the cross-sections that need to be flipped, according to the
        % flip_sections variable
        FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
        flip_notches(FlipName,AllSectionsName,FlippedOutputName);

    else
        disp('A flip index vector for the chosen data has been found.');
    end
    
    load(FlippedOutputName);
    
    

    %% Get ellipse fits and difference data 
    % Check to see if there is already an ellipse fit for the chosen distance
    
    % Get file name to look for (same basic code as before, but looking for
    % existing ellipse fit data instead of flipping data)
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

    if ~isfile(searchname)
        disp('No ellipse data for this slice distance exists in the current folder. Create data now.');

        % Calculate the ellipse fits for all cross-sections at chosen distance
        AllEllipseName = strcat(output_prefix,'_AllEllipses.mat');
        ellipse_fitting_V2(FlippedOutputName,AllEllipseName);

        load(AllEllipseName);
        load(FlippedOutputName);

        % Plot the ellipse fits and see if any of them have problems. Tag
        % the problem cross-sections with a 1, as with the flip data.
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
            else
                continue
            end    
        end    
        
        % This section gives the user a chance to correct any issues with
        % the ellipse fits. Sometimes the ellipse fitting will be very
        % noisy or an ellipse won't be fit, showing a huge spiral. These
        % cases MUST be rejected for PCA to work properly.
        while 1
            fixes_needed = input('Does the ellipse problems vector need manual correction? Y/N ','s');
            switch fixes_needed
                case 'Y'
                    openvar('problem_indices');
                    disp('Giving control to keyboard for manual editing of ellipse data.');
                    disp('Use dbcont command to exit keyboard editing mode.');
                    keyboard;
                    break
                case 'N'
                    break

                otherwise
                    disp('Not a recognized answer. Please try again.');
            end
        end

        % Save problem_indices for later use
        ProblemEllipses = strcat(output_prefix,'_ProblemEllipses.mat');
        save(ProblemEllipses,'problem_indices');

        problem_indices

    % This code runs if there is existing ellipse data in the working
    % folder. This prevents the user from having to do the tedious sorting
    % tasks every time they want to run the data.
    else
        disp('Ellipse data for the chosen slice distance has been found.');
        AllEllipseName = strcat(output_prefix,'_AllEllipses.mat');

        try
            ProblemEllipses = strcat(output_prefix,'_ProblemEllipses.mat');
            load(ProblemEllipses,'problem_indices');
        catch
            ProblemEllipses = strcat('Random_',output_prefix,'_ProblemEllipses.mat');
            load(ProblemEllipses,'problem_indices');
        end
    end
    
    
    %% Prepare data for PCA
    % Add the difference data and corresponding ellipse fits to separate
    % large arrays for feeding into PCA. Take note of the starting and
    % ending indices corresponding to each slice distance so data can be
    % correctly reconstructed later. Also save the error_indices that
    % correspond.
    try
        AllEllipseName = strcat(output_prefix,'_AllEllipses.mat');
        load(AllEllipseName);
        load(ProblemEllipses,'problem_indices');
    catch
        AllEllipseName = strcat('Random_',output_prefix,'_AllEllipses.mat');
        load(AllEllipseName);
        load(ProblemEllipses,'problem_indices');
    end
    
    % Gathering some size parameters for later use
    n_ALL_DIFF_R_ext      = size(ALL_DIFF_R_ext,1);
    n_ALL_DIFF_R_int      = size(ALL_DIFF_R_int,1);
    n_ALL_R_ext           = size(ALL_R_ext,1);
    n_ALL_R_int           = size(ALL_R_int,1);
    n_ALL_A               = size(ALL_A,1);
    n_ALL_B               = size(ALL_B,1);
    n_ALL_ELLIPSE_T       = size(ALL_ELLIPSE_T,1);
    n_ALL_AVG_RIND_T      = size(ALL_AVG_RIND_T,1);
    n_ALL_ELLIPSE_R_ext   = size(ALL_ELLIPSE_R_ext);
    n_ALL_ELLIPSE_R_int   = size(ALL_ELLIPSE_R_int);
    
    % This index is to determine the row for arrays such as
    % slice_startstop, ALL_PROBLEM_INDICES, etc.
    ind = find(slice_dists == slice);
    
    ALL_PROBLEM_INDICES(ind,:) = {problem_indices}; % EACH ROW IS A NEW SLICE DISTANCE
    
        
    % Remove error cases and adjust indices, starting at the bottom of the
    % array. Save adjusted indices.
    % There are normally 980 slices in the data set that I used. Initialize
    % this variable with all the expected slice numbers.
    adj_indices_temp = linspace(1,980,980);
    
    % Get rid of any indices that had a problem.
    adj_indices_temp(problem_indices) = []; % to map between stalk number and index value, use the index and the value of adj_indices
    
    % Save this shifted collection of indices for reconstructing specific
    % cross-sections later on (in transverse_wrapper_V4, for example)
    adj_indices{ind,1} = adj_indices_temp;
    
    % Remove data at problem_indices for the current slice location
    DIFF_R_ext(problem_indices,:) = [];
    DIFF_R_int(problem_indices,:) = [];
    R_ext(problem_indices,:) = [];
    R_int(problem_indices,:) = [];
    A(problem_indices,:) = [];
    B(problem_indices,:) = [];
    ELLIPSE_T(problem_indices,:) = [];
    AVG_RIND_T(problem_indices,:) = [];
    ELLIPSE_R_ext(problem_indices,:) = [];
    ELLIPSE_R_int(problem_indices,:) = [];
    
    % Add the data from the current slice to the "ALL" version of each data
    % array for saving at the end of this script and sending to PCA.
    ALL_DIFF_R_ext      = [ALL_DIFF_R_ext; DIFF_R_ext];
    ALL_DIFF_R_int      = [ALL_DIFF_R_int; DIFF_R_int];
    ALL_R_ext           = [ALL_R_ext; R_ext];
    ALL_R_int           = [ALL_R_int; R_int];
    ALL_A               = [ALL_A; A];
    ALL_B               = [ALL_B; B];
    ALL_ELLIPSE_T       = [ALL_ELLIPSE_T; ELLIPSE_T];
    ALL_AVG_RIND_T      = [ALL_AVG_RIND_T; AVG_RIND_T];
    ALL_ELLIPSE_R_ext   = [ALL_ELLIPSE_R_ext; ELLIPSE_R_ext];
    ALL_ELLIPSE_R_int   = [ALL_ELLIPSE_R_int; ELLIPSE_R_int];
    
    % Since all the data is being lumped together from lots of different
    % slice locations, we need a way to later determine which slice
    % location a given cross-section comes from. slice_startstop gives the
    % slice location and the starting and ending indices for that slice
    % location in the "ALL" data variables. This will be useful later.
    if ind == 1
        slice_startstop(ind,2) = 1;
    else
        slice_startstop(ind,2) = n_ALL_DIFF_R_ext + 1;
    end
    
    slice_startstop(ind,3) = size(ALL_DIFF_R_ext,1);
        
    
    %% Check the individual slice distance PCA results to make sure nothing was missed
    % Run PCA on the resulting large data set. Save this for access by
    % transverse_wrapper_V4.m
    [ext_rhoPCAs, ext_rhocoeffs, ext_rhoPCA_variances, ext_rhotstat, ext_rhoexplained, ext_rhovarMeans] = pca(ALL_DIFF_R_ext,'Centered',false);
    [int_rhoPCAs, int_rhocoeffs, int_rhoPCA_variances, int_rhotstat, int_rhoexplained, int_rhovarMeans] = pca(ALL_DIFF_R_int,'Centered',false);
    
    % Create arrays for "explained" data, which shows the successive
    % contributions each principal component makes to the approximation.
    % These data will be used to create the plots below.
    ext_rhoexplained_tot = zeros(size(ext_rhoexplained));
    int_rhoexplained_tot = zeros(size(ext_rhoexplained));
    for i = 1:length(ext_rhoexplained_tot)
        ext_rhoexplained_tot(i) = sum(ext_rhoexplained(1:i));
        int_rhoexplained_tot(i) = sum(int_rhoexplained(1:i));
    end
    
    % Plot the exterior explained data (current slice location only)
    figure(1);
    plot(ext_rhoexplained_tot(1:50,1),'-*');
    title('Exterior Non-Ellipse PCs');
    xlabel('# of PCs');
    ylabel('% Variance Explained');
    
    % Plot the interior explained data (current slice location only)
    figure(2);
    plot(int_rhoexplained_tot(1:50,1),'-*');
    title('Interior Non-Ellipse PCs');
    xlabel('# of PCs');
    ylabel('% Variance Explained');
    
    ELLIPSE_T = ELLIPSE_T';
    theta = ELLIPSE_T(:,1);
    ELLIPSE_T = ELLIPSE_T';
    
    % Polar plot of the exterior principal components (current slice location only)
    figure(3);
    polarplot(theta,ext_rhoPCAs(:,1));
    hold on
    polarplot(theta,ext_rhoPCAs(:,2));
    polarplot(theta,ext_rhoPCAs(:,3));
    polarplot(theta,ext_rhoPCAs(:,4));
    polarplot(theta,ext_rhoPCAs(:,5));
    title('Exterior Rho Principal Components');
    legend('PC1','PC2','PC3','PC4','PC5');
    
    % Polar plot of the interior principal components (current slice location only)
    figure(4);
    polarplot(theta,int_rhoPCAs(:,1));
    hold on
    polarplot(theta,int_rhoPCAs(:,2));
    polarplot(theta,int_rhoPCAs(:,3));
    polarplot(theta,int_rhoPCAs(:,4));
    polarplot(theta,int_rhoPCAs(:,5));
    title('Interior Rho Principal Components');
    legend('PC1','PC2','PC3','PC4','PC5');

    pause();
    
    
    
    
end

%% Run PCA on combined data set
% Run PCA on the resulting large data set. Save this for access by
% transverse_wrapper_V4.m
[ext_rhoPCAs, ext_rhocoeffs, ext_rhoPCA_variances, ext_rhotstat, ext_rhoexplained, ext_rhovarMeans] = pca(ALL_DIFF_R_ext,'Centered',false);
[int_rhoPCAs, int_rhocoeffs, int_rhoPCA_variances, int_rhotstat, int_rhoexplained, int_rhovarMeans] = pca(ALL_DIFF_R_int,'Centered',false);

% Calculate the variance explained progression for plotting below
ext_rhoexplained_tot = zeros(size(ext_rhoexplained));
int_rhoexplained_tot = zeros(size(ext_rhoexplained));
for i = 1:length(ext_rhoexplained_tot)
    ext_rhoexplained_tot(i) = sum(ext_rhoexplained(1:i));
    int_rhoexplained_tot(i) = sum(int_rhoexplained(1:i));
end

% Plot the exterior explained data (all slice locations)
figure(1);
plot(ext_rhoexplained_tot(1:50,1),'-*');
title('ALL DATA Exterior Non-Ellipse PCs');
xlabel('# of PCs');
ylabel('% Variance Explained');

% Plot the interior explained data (all slice locations)
figure(2);
plot(int_rhoexplained_tot(1:50,1),'-*');
title('ALL DATA Interior Non-Ellipse PCs');
xlabel('# of PCs');
ylabel('% Variance Explained');

ELLIPSE_T = ELLIPSE_T';
theta = ELLIPSE_T(:,1);
ELLIPSE_T = ELLIPSE_T';

% Polar plot of the exterior principal components (all slice locations)
figure(3);
polarplot(theta,ext_rhoPCAs(:,1));
hold on
polarplot(theta,ext_rhoPCAs(:,2));
polarplot(theta,ext_rhoPCAs(:,3));
polarplot(theta,ext_rhoPCAs(:,4));
polarplot(theta,ext_rhoPCAs(:,5));
title('ALL DATA Exterior Rho Principal Components');
legend('PC1','PC2','PC3','PC4','PC5');

% Polar plot of the interior principal components (all slice locations)
figure(4);
polarplot(theta,int_rhoPCAs(:,1));
hold on
polarplot(theta,int_rhoPCAs(:,2));
polarplot(theta,int_rhoPCAs(:,3));
polarplot(theta,int_rhoPCAs(:,4));
polarplot(theta,int_rhoPCAs(:,5));
title('ALL DATA Interior Rho Principal Components');
legend('PC1','PC2','PC3','PC4','PC5');




%% Save the final data in a new mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, SaveName);
save(SaveFile,'ELLIPSE_T','ELLIPSE_R_ext','ext_rhocoeffs',...
    'ext_rhoPCAs','ext_rhoexplained','ext_rhovarMeans','ELLIPSE_R_int',...
    'int_rhocoeffs','int_rhoPCAs','int_rhoexplained','int_rhovarMeans',...
    'ALL_DIFF_R_ext','ALL_DIFF_R_int','ALL_ELLIPSE_T','ALL_ELLIPSE_R_ext',...
    'ALL_ELLIPSE_R_int','slice_startstop','ALL_PROBLEM_INDICES',...
    'adj_indices','slice_dists','ALL_A','ALL_B','ALL_R_ext','ALL_R_int',...
    'ALL_AVG_RIND_T');


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Localized functions
function ChooseSections(method,stalknums,dist,Table,error_indices,npoints,SaveName)
% FILENAME: ChooseSections.m
% AUTHOR: Ryan Larson
% DATE: 7/3/2019
%
% PURPOSE: The slice locations were created by the CT scan process, so they
% aren't all at "nice" distances. This function selects the desired
% cross-sections based on a method. 
% 
% 
% INPUTS:
%       method: String that determines how the cross-sections will be
%       selected. There are two methods: 'samedist' and 'all'. However, the
%       way later code was constructed means the only method that was used
%       was 'samedist'.
% 
%       stalknums: A row vector of integers corresponding to the stalk
%       numbers. These integers will be somewhere between 1 and 980.
% 
%       dist: A number (integer or float) that determines the desired slice
%       location, relative to the node. It can be positive or negative, and
%       is measured in mm. Maximum values should be approximately +/- 40mm.
% 
%       Table: A table containing the downsampled, centered, and rotated
%       cross-section boundaries. This is the table that is output from
%       Jared's script. Usually this will be Stalk_TableDCR, which loads as
%       a variable in StalksDCR_360pts_V2.mat.
% 
%       error_indices: Indices that had problems in pre-processing and
%       should be ignored. Also loads as a variable in
%       StalksDCR_360pts_V2.mat.
% 
%       npoints: The number of sample points the cross-sections have been
%       downsampled to. This was previously determined by Jared's
%       preprocessing script, and is only here to carry over the proper
%       indexing.
% 
%       SaveName: String name to save the data to as a .mat file. Must
%       include .mat in the name.
%       
% OUTPUTS:
%       N/A
%
%
% NOTES: 
%       Other outputs:
%           - Named .mat file that contains all the ellipse fit and PCA
%           data, which will be used when running transverse_wrapper_V4.m.
%       
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

allrows = size(Table,1);

switch method
    
    % Choose a number of cross-sections that are all at the same distance
    % from the node
    case 'samedist'
        indices = zeros(1,length(stalknums));

        
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
            'int_T','int_Rho','avg_rind_thick','indices','selectedTable',...
            'npoints','deletesections');

        
    % Choose all cross-sections (not really useful in hindsight, but was
    % written before full PCA process was defined)
    case 'all'        
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
    
    % Mark any cross-sections that are made up of zeros by setting them to
    % something that will be obviously wrong when manually checking ellipse
    % fits
    
    if (all(xcut == 0) || all(ycut == 0))
        xcut = 1000*[11.9371843338013;11.8635606765747;11.8245401382446;11.7747907638550;11.7636137008667;11.7513217926025;11.7293729782105;11.7618465423584;11.7397794723511;11.8146848678589;11.7211742401123;11.6160907745361;11.5834693908691;11.5488634109497;11.5280466079712;11.4122400283813;11.2877616882324;11.2347469329834;11.1380500793457;10.9807081222534;10.8686056137085;10.7776327133179;10.6507196426392;10.5666399002075;10.4331426620483;10.2808752059937;10.1859483718872;10.0852527618408;9.95389461517334;9.82945919036865;9.72529697418213;9.61997604370117;9.46991920471191;9.31812000274658;9.16870021820068;9.04158782958984;8.90647315979004;8.77705669403076;8.67132759094238;8.53238964080811;8.35374259948731;8.16121101379395;7.98431634902954;7.85418319702148;7.72590875625610;7.58815336227417;7.45849895477295;7.30975627899170;7.13432264328003;6.97212791442871;6.83072328567505;6.67568683624268;6.50259399414063;6.34706783294678;6.19373226165772;6.02386379241943;5.83662033081055;5.66628932952881;5.51614332199097;5.36890840530396;5.19331264495850;5.03751039505005;4.87642765045166;4.69888496398926;4.51799726486206;4.34112358093262;4.18711233139038;4.01929712295532;3.82844710350037;3.65820479393005;3.49178791046143;3.31141233444214;3.14243388175964;2.98499917984009;2.80316996574402;2.63344144821167;2.45546531677246;2.27513217926025;2.09673070907593;1.92558479309082;1.75083839893341;1.57933926582336;1.40231263637543;1.22391927242279;1.04622578620911;0.871519982814789;0.697965204715729;0.523268401622772;0.348291426897049;0.174448028206825;-4.37227299698861e-07;-0.174284666776657;-0.348030894994736;-0.520720541477203;-0.693446636199951;-0.865429878234863;-1.03545892238617;-1.20444118976593;-1.37653625011444;-1.54755997657776;-1.71807074546814;-1.89003670215607;-2.05707454681397;-2.21517491340637;-2.37523984909058;-2.55128669738770;-2.72201323509216;-2.89396500587463;-3.07621908187866;-3.24406313896179;-3.40437126159668;-3.55446076393127;-3.70386385917664;-3.85397577285767;-4.01397371292114;-4.18080663681030;-4.34828853607178;-4.50251483917236;-4.65524196624756;-4.84495019912720;-4.99901247024536;-5.14868736267090;-5.31385993957520;-5.45570230484009;-5.62545299530029;-5.78521919250488;-5.93097782135010;-6.07301187515259;-6.20860242843628;-6.35488939285278;-6.55632448196411;-6.72821474075317;-6.89345264434814;-7.08150720596314;-7.24557924270630;-7.35573911666870;-7.17333126068115;-6.94547891616821;-6.79845809936523;-6.61049652099609;-6.43506336212158;-6.27632665634155;-6.13361406326294;-5.94674444198608;-5.77380752563477;-5.62645292282105;-5.48189973831177;-5.29240417480469;-5.12724733352661;-4.99361371994019;-4.84064579010010;-4.64613723754883;-4.47754716873169;-4.31419754028320;-4.14768457412720;-3.97254443168640;-3.82159447669983;-3.66316986083984;-3.48787474632263;-3.30496263504028;-3.15222120285034;-2.99911570549011;-2.83224821090698;-2.66091227531433;-2.50481224060059;-2.34091997146606;-2.16448307037354;-1.99323868751526;-1.82688367366791;-1.65910267829895;-1.49770700931549;-1.33571767807007;-1.16655278205872;-1.00036275386810;-0.837868511676788;-0.669154405593872;-0.500231862068176;-0.333453625440598;-0.166822329163551;1.13528372480687e-07;0.166331484913826;0.334067136049271;0.501357376575470;0.666979789733887;0.832698047161102;0.998343646526337;1.16228866577148;1.32627189159393;1.49585855007172;1.65915548801422;1.82434475421906;1.99517750740051;2.16472029685974;2.33758425712585;2.50638127326965;2.68204832077026;2.85898113250732;3.02217745780945;3.18090009689331;3.34416341781616;3.51388120651245;3.69013237953186;3.85169339179993;4.00911951065064;4.20240354537964;4.36277866363525;4.52431869506836;4.68642330169678;4.84755754470825;5.03232288360596;5.20052623748779;5.37546968460083;5.56187343597412;5.72971630096436;5.87029695510864;6.03825664520264;6.18831110000610;6.33646345138550;6.51382923126221;6.69243431091309;6.84649372100830;7.01449489593506;7.16981601715088;7.31131458282471;7.47539043426514;7.63156509399414;7.85899925231934;8.01014423370361;8.19667720794678;8.35878944396973;8.47652244567871;8.60494613647461;8.75809860229492;8.94601535797119;9.08240222930908;9.20913982391357;9.40316677093506;9.52799606323242;9.63124275207520;9.80542659759522;9.93904399871826;10.0585393905640;10.1908664703369;10.3225965499878;10.4634752273560;10.5509290695190;10.6960325241089;10.7989549636841;10.8371763229370;10.9699687957764;11.0869083404541;11.2064266204834;11.3318710327148;11.3706150054932;11.4170074462891;11.5129404067993;11.5899295806885;11.6690845489502;11.7415847778320;11.7823114395142;11.7643365859985;11.8403053283691;11.8786678314209;11.9376668930054;11.9972867965698;11.9749670028687;12.0239715576172;11.9950208663940;11.9469957351685];
        ycut = 100*[0;0.207079216837883;0.412922024726868;0.617090642452240;0.822591960430145;1.02810728549957;1.23280680179596;1.44417321681976;1.64991843700409;1.87126231193542;2.06675934791565;2.25793933868408;2.46214246749878;2.66626501083374;2.87426495552063;3.05790042877197;3.23671364784241;3.43480658531189;3.61897182464600;3.78096103668213;3.95584869384766;4.13714599609375;4.30317020416260;4.48527240753174;4.64513444900513;4.79405117034912;4.96801900863648;5.13869333267212;5.29258012771606;5.44855785369873;5.61490297317505;5.78026485443115;5.91746234893799;6.05125808715820;6.18436574935913;6.33098793029785;6.47093152999878;6.61398601531982;6.77478361129761;6.90939235687256;7.00962162017822;7.09443235397339;7.18911075592041;7.32414388656616;7.46082401275635;7.58815336227417;7.72350168228149;7.83875370025635;7.92346906661987;8.02051544189453;8.14053916931152;8.24378681182861;8.32293987274170;8.42284393310547;8.52494144439697;8.60296916961670;8.65314483642578;8.72532176971436;8.82767391204834;8.93536472320557;8.99508285522461;9.08791065216065;9.17122554779053;9.22208023071289;9.26326656341553;9.30957126617432;9.40440940856934;9.46887207031250;9.47573757171631;9.52994823455811;9.59360790252686;9.61703968048096;9.67141819000244;9.76349449157715;9.77581405639648;9.82813549041748;9.84833240509033;9.85468101501465;9.86434459686279;9.90627765655518;9.92949485778809;9.97155475616455;9.97797203063965;9.96802330017090;9.95417690277100;9.96151161193848;9.98136329650879;9.98455238342285;9.97376155853272;9.99413394927979;10.0025949478149;9.98479270935059;9.96631050109863;9.93594264984131;9.91674804687500;9.89190578460693;9.85172939300537;9.80939102172852;9.79456615447998;9.77091026306152;9.74366283416748;9.72339439392090;9.67777252197266;9.59497928619385;9.52656745910645;9.52153205871582;9.49278831481934;9.46573066711426;9.46762752532959;9.42144489288330;9.35343360900879;9.25968742370606;9.16738414764404;9.07939815521240;9.01553153991699;8.96576690673828;8.91531276702881;8.83668041229248;8.75523567199707;8.74052238464356;8.65854263305664;8.56885623931885;8.50395202636719;8.40104484558106;8.34007835388184;8.26214885711670;8.16329193115234;8.05915737152100;7.94664907455444;7.84763669967651;7.81352329254150;7.73992681503296;7.65595388412476;7.59398746490479;7.50301551818848;-7.61709213256836;-7.69245624542236;-7.71373748779297;-7.82073163986206;-7.87808513641357;-7.94664192199707;-8.03333091735840;-8.13957786560059;-8.18499565124512;-8.24585437774658;-8.34156036376953;-8.44138431549072;-8.46961498260498;-8.53317642211914;-8.64919471740723;-8.73275661468506;-8.73811149597168;-8.78767681121826;-8.84541988372803;-8.89474010467529;-8.92248153686523;-9.00310993194580;-9.06665802001953;-9.08623027801514;-9.08031272888184;-9.15471458435059;-9.23032569885254;-9.26385974884033;-9.27970981597900;-9.34808826446533;-9.38891601562500;-9.37540149688721;-9.37744045257568;-9.39850902557373;-9.40924167633057;-9.45614624023438;-9.50411415100098;-9.50082778930664;-9.51782798767090;-9.57688426971436;-9.56934356689453;-9.54496097564697;-9.54891967773438;-9.55730628967285;-9.52029418945313;-9.52917194366455;-9.56648254394531;-9.56643390655518;-9.53824234008789;-9.51778316497803;-9.49861526489258;-9.46609783172607;-9.43690299987793;-9.44447422027588;-9.40953922271729;-9.38544654846191;-9.38658237457275;-9.37642669677734;-9.37553691864014;-9.35394382476807;-9.35341930389404;-9.35131454467773;-9.30130100250244;-9.23800373077393;-9.18801498413086;-9.15397739410400;-9.13339233398438;-9.07401752471924;-9.00462818145752;-9.01208591461182;-8.94502544403076;-8.87947082519531;-8.81387805938721;-8.74522590637207;-8.71624088287354;-8.65513324737549;-8.60254573822022;-8.56453132629395;-8.49465370178223;-8.38365554809570;-8.31095123291016;-8.21216297149658;-8.11030197143555;-8.04391002655029;-7.97573423385620;-7.87599372863770;-7.79038429260254;-7.68868541717529;-7.57108879089356;-7.47539281845093;-7.36971330642700;-7.32863378524780;-7.21236515045166;-7.12526369094849;-7.01386022567749;-6.86414957046509;-6.72291898727417;-6.59969997406006;-6.49966144561768;-6.35956859588623;-6.21164035797119;-6.10648632049561;-5.95375299453735;-5.78703594207764;-5.66116809844971;-5.50929927825928;-5.34821939468384;-5.19250631332398;-5.03466796875000;-4.87920141220093;-4.69757413864136;-4.54019594192505;-4.36306142807007;-4.16000366210938;-3.99274492263794;-3.81752705574036;-3.64118814468384;-3.46450138092041;-3.26047325134277;-3.05917525291443;-2.87049674987793;-2.67574572563171;-2.48034143447876;-2.28233480453491;-2.07753682136536;-1.86328649520874;-1.66404628753662;-1.45851802825928;-1.25470173358917;-1.04962432384491;-0.837370157241821;-0.630149722099304;-0.418876588344574;-0.208537995815277];
    end
    
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
