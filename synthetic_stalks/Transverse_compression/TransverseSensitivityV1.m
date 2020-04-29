function TransverseSensitivityV1(slices,stalknums,AllSlicesPCA,percent_change,numNEPCs,plotting)
% FILENAME: TransverseSensitivityV1.m
% AUTHOR: Ryan Larson
% DATE: 1/24/2020
%
% PURPOSE: Create Python scripts to be fed into Abaqus for a sensitivity
% study. This script can be used for any of the cross-sections in the PCA
% data, depending on values chosen for slices and stalknums.
% 
% 
% INPUTS:
%       slices:  A subset of the input that went into AllTransversePCA.m
%       (as of 1/22/2020, this was
%       [-40 -30 -20 -15 -10 -5 0 5 10 15 20 30 40])
% 
%       stalknums: A vector of unique integers from 1 to 980 that
%       determines which stalks to sample from (use randperm(980,K) to
%       choose K unique integers from 1 to 980)
% 
%       AllSlicesPCA: PCA data output from AllTransversePCa.m. Enter this
%       as a string ('AllSlicesPCA.mat').
% 
%       percent_change: The percentage to change each parameter value by
%       for the sensitivity study. Enter as a decimal (i.e. 10% would be
%       0.1)
% 
%       numNEPCs: Number of principal components to include in the base
%       model (Ryan used 5)
% 
%       plotting: Enter 1 to plot, 0 to not plot. This is useful for
%       checking outputs.
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
%       - 
% 
% -------------------------------------------------------------------------
% SUBROUTINES:
%   make_case.m: Create a specimen-specific, model-specific Python
%   script that runs the chosen model in transverse compression when
%   fed to ABAQUS.
% 
%   writespline_V2.m: A subroutine of make_case.m. Convert boundary
%   data from Matlab arrays to strings that can be used in the Python
%   scripts.
% 
%   rpts.m: Create an ellipse from major and minor diameter values.
% 
%   get_materials.m: Generate rind and pith stiffnesses to use in a
%   given model.
% 
% PSEUDO-CODE:
%   Load PCA data.
%   Determine factor for multiplying parameters by in sensitivity study.
%   
%   for each slice distance:
%       Find the place in the PCA and ellipse data where the data for the
%       current slice starts.
%       
%       for each stalk:
%           Get the actual index for the chosen slice (adjusting for cases
%           that were deleted due to bad ellipse fits or otherwise).
% 
%           Generate a base approximation profile (base case) using ellipse
%           fit data the first five principal components.
% 
%           Generate rind and pith material properties to be used for each
%           sensitivity case.
% 
%           Create a Python script for the base case (case 0).
% 
%           Calculate a new ellipse fit with the major diameter stretched
%           by the sensitivity factor (all other parameters the same as the
%           base case).
%           Create a Python script for the major diameter sensitivity case
%           (case 1).
% 
%           Calculate a new ellipse fit with the minor diameter stretched
%           by the sensitivity factor (all other parameters the same as the
%           base case).
%           Create a Python script for the minor diameter sensitivity case
%           (case 2).
% 
%           Use the base case exterior profile to calculate a new ellipse
%           interior using the base rind thickness multiplied by the
%           sensitivity factor (all other parameters the same as the base
%           case).
%           Create a Python script for the rind thickness sensitivity case
%           (case 3).
% 
%           Use the base case geometry, and multiply the rind stiffness by
%           the sensitivity factor (all other parameters the same as the
%           base case).
%           Create a Python script for the rind stiffness sensitivity case
%           (case 4).
% 
%           Use the base case geometry, and multiply the pith stiffness by
%           the sensitivity factor (all other parameters the same as the
%           base case).
%           Create a Python script for the pith stiffness sensitivity case
%           (case 5).
%   
%           Use the base case ellipse parameters and add in the principal
%           component features, with the 1st principal component scaling
%           factor multiplied by the sensitivity factor (all other
%           parameters the same as the base case).
%           Create a Python script for the PC 1 sensitivity case (case 6).
% 
%           Use the base case ellipse parameters and add in the principal
%           component features, with the 3nd principal component scaling
%           factor multiplied by the sensitivity factor (all other
%           parameters the same as the base case).
%           Create a Python script for the PC 2 sensitivity case (case 7).
% 
%           Use the base case ellipse parameters and add in the principal
%           component features, with the 3rd principal component scaling
%           factor multiplied by the sensitivity factor (all other
%           parameters the same as the base case).
%           Create a Python script for the PC 3 sensitivity case (case 8).
% 
%           Use the base case ellipse parameters and add in the principal
%           component features, with the 4th principal component scaling
%           factor multiplied by the sensitivity factor (all other
%           parameters the same as the base case).
%           Create a Python script for the PC 4 sensitivity case (case 9).
% 
%           Use the base case ellipse parameters and add in the principal
%           component features, with the 5th principal component scaling
%           factor multiplied by the sensitivity factor (all other
%           parameters the same as the base case).
%           Create a Python script for the PC 5 sensitivity case (case 10).
% 
%       end
% 
%   end
%           
% 
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - General sensitivity study, best applied to all data being examined
% with transverse_wrapper_V4.m
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

%% Initial variables
set(0,'DefaultFigureWindowStyle','docked');
load(AllSlicesPCA);
group = 1;
problem_slice_stalk = [];

% Get multiplier based on the desired percent change for each parameter
plus_change = 1 + percent_change;

write_Python_template3;  % Create Template cell array that can be copied and used to make individualized Python scripts


%% Create all geometry cases for a given cross section
% Iterate through slices (determine group number here)
for slice = slices
    
    % Determine the row in the PCA data where the data for the current
    % slice starts.
    sliceidx = find(slice_dists == slice);
    startidx = slice_startstop(sliceidx,2);
    
    % For each slice position, iterate through the stalks
    for stalk = stalknums
        % Get the actual index of the chosen data and create a Python script for
        % that case, numbering by group
        indices = cell2mat(adj_indices(sliceidx,1));
        stalkidx = find(indices == stalk);
        
        % If the chosen stalk had a bad ellipse fit, catch this so further
        % cases can be generated if needed.
        if isempty(stalkidx)
            problem_slice_stalk = [problem_slice_stalk; slice, stalk];
            continue
        end
        
        % Calculate the row index corresponding to the current specific
        % cross-section, as found in the PCA data
        adj_ind = startidx + stalkidx - 1;
        
        %% Create case from ellipse and PCA data (using "ALL" variables)
        GROUP = sprintf('%d',group); % Group number
        ID = sprintf('%d',stalk); % Cross-section number
        
        write_Python_template3;
        
        % Get random material properties
        material_method = 'random';
        [Erind,Epith] = get_materials('random');
        
        % Make profile that includes the chosen number of principal
        % components
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:numNEPCs
            % Add all PCs (through the current PC) to the ellipse in polar coordinates
            NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
        end
        
        % Calculate the base case exterior and interior profiles
        base_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
        base_int = normintV2(base_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));
        
        %% Base case (case 0)
        % Write the Python script for the base case
        case_num = 0;
        Script = Template; % Reset the script template
        make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,ALL_ELLIPSE_T,Script,Erind,Epith)
        
        %% Change A (case 1)
        % Write the Python script for case 1
        % Adjust exterior points in polar
        Tnew = ALL_ELLIPSE_T(adj_ind,:);
        
        % Stretch the major diameter of the ellipse fit by percent_change
        Aplus_ellipse = rpts(360,ALL_ELLIPSE_T(adj_ind,:),(plus_change*ALL_A(adj_ind)),ALL_B(adj_ind));
        
        % Create the new profile using the new ellipse
        Aplus_ext = Aplus_ellipse - NEPC_ext;

        % Calculate the interior points
        Aplus_int = normintV2(Aplus_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

        % Check shape 
        if plotting == 1
            polarplot(Tnew,Aplus_ext,'r');
            hold on
            polarplot(Tnew,Aplus_int,'r');
            polarplot(Tnew,base_ext,'b');
            polarplot(Tnew,base_int,'b');
            title('Changing A');
            pause();
            close;
        end
        
        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,adj_ind,ID,GROUP,Aplus_ext,Aplus_int,Tnew,Script,Erind,Epith);
        

        %% Change B (case 2)
        % Write the Python script for case 2
        % Adjust exterior points in polar
        Tnew = ALL_ELLIPSE_T(adj_ind,:);
        
        % Stretch the minor diameter of the ellipse fit by percent_change
        Bplus_ellipse = rpts(360,ALL_ELLIPSE_T(adj_ind,:),(ALL_A(adj_ind)),plus_change*ALL_B(adj_ind));
        
        % Create the new profile using the new ellipse
        Bplus_ext = Bplus_ellipse - NEPC_ext;

        % Calculate the interior points
        Bplus_int = normintV2(Bplus_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

        % Check shape
        if plotting == 1
            polarplot(Tnew,Bplus_ext,'r');
            hold on
            polarplot(Tnew,Bplus_int,'r');
            polarplot(Tnew,base_ext,'b');
            polarplot(Tnew,base_int,'b');
            title('Changing B');
            pause();
            close;
        end

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,adj_ind,ID,GROUP,Bplus_ext,Bplus_int,Tnew,Script,Erind,Epith);

        
        %% Change T (case 3)
        Tnew = ALL_ELLIPSE_T(adj_ind,:);

        % Calculate the interior points
        Tplus_int = normintV2(base_ext,ALL_ELLIPSE_T(adj_ind,:),plus_change*ALL_AVG_RIND_T(adj_ind));

        % Check shape
        if plotting == 1
            polarplot(Tnew,base_ext,'r');
            hold on
            polarplot(Tnew,Tplus_int,'r');
            polarplot(Tnew,base_ext,'b');
            polarplot(Tnew,base_int,'b');
            title('Changing T');
            pause();
            close;
        end

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,adj_ind,ID,GROUP,base_ext,Tplus_int,Tnew,Script,Erind,Epith);        
        
        %% Change Erind (case 4)
        Tnew = ALL_ELLIPSE_T(adj_ind,:);
        
        % Calculate the new Erind
        Erind_plus = Erind*plus_change;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,Tnew,Script,Erind_plus,Epith);        
        
        %% Change Epith (case 5)
        Tnew = ALL_ELLIPSE_T(adj_ind,:);
        
        % Calculate the new Erind
        Epith_plus = Epith*plus_change;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,Tnew,Script,Erind,Epith_plus);        
        
        %% Change NEPC 1 (case 6)
        if numNEPCs >= 1
        
            Tnew = ALL_ELLIPSE_T(adj_ind,:);

            % Make profile that includes NEPCs 1-5
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

            for k = 1:numNEPCs
                if k == 1
                    NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                else
                    NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                end
            end

            R_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            R_int = normintV2(R_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

            % Check shape
            if plotting == 1
                polarplot(Tnew,R_ext,'r');
                hold on
                polarplot(Tnew,R_int,'r');
                polarplot(Tnew,base_ext,'b');
                polarplot(Tnew,base_int,'b');
                title('Changing NEPC 1');
                pause();
                close;
            end

            % Create cases
            case_num = case_num + 1;
            Script = Template; % Reset the script template    
            make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,Tnew,Script,Erind,Epith);        
            
        end
        
        %% Change NEPC 2 (case 7)
        if numNEPCs >= 2
            
            Tnew = ALL_ELLIPSE_T(adj_ind,:);

            % Make profile that includes NEPCs 1-5
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

            for k = 1:numNEPCs
                if k == 2
                    NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                else
                    NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                end
            end

            R_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            R_int = normintV2(R_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

            % Check shape
            if plotting == 1
                polarplot(Tnew,R_ext,'r');
                hold on
                polarplot(Tnew,R_int,'r');
                polarplot(Tnew,base_ext,'b');
                polarplot(Tnew,base_int,'b');
                title('Changing NEPC 2');
                pause();
                close;
            end

            % Create cases
            case_num = case_num + 1;
            Script = Template; % Reset the script template    
            make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,Tnew,Script,Erind,Epith);        
        
        end
        
        %% Change NEPC 3 (case 8)
        if numNEPCs >= 3
            
            Tnew = ALL_ELLIPSE_T(adj_ind,:);

            % Make profile that includes NEPCs 1-5
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

            for k = 1:numNEPCs
                if k == 3
                    NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                else
                    NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                end
            end

            R_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            R_int = normintV2(R_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

            % Check shape
            if plotting == 1
                polarplot(Tnew,R_ext,'r');
                hold on
                polarplot(Tnew,R_int,'r');
                polarplot(Tnew,base_ext,'b');
                polarplot(Tnew,base_int,'b');
                title('Changing NEPC 3');
                pause();
                close;
            end

            % Create cases
            case_num = case_num + 1;
            Script = Template; % Reset the script template    
            make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,Tnew,Script,Erind,Epith);        
        
        end
        
        %% Change NEPC 4 (case 9)
        if numNEPCs >= 4
            
            Tnew = ALL_ELLIPSE_T(adj_ind,:);

            % Make profile that includes NEPCs 1-5
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

            for k = 1:numNEPCs
                if k == 4
                    NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                else
                    NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                end
            end

            R_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            R_int = normintV2(R_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

            % Check shape
            if plotting == 1
                polarplot(Tnew,R_ext,'r');
                hold on
                polarplot(Tnew,R_int,'r');
                polarplot(Tnew,base_ext,'b');
                polarplot(Tnew,base_int,'b');
                title('Changing NEPC 4');
                pause();
                close;
            end

            % Create cases
            case_num = case_num + 1;
            Script = Template; % Reset the script template    
            make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,Tnew,Script,Erind,Epith);        
        
        end
        
        %% Change NEPC 5 (case 10)
        if numNEPCs >= 5
            
            Tnew = ALL_ELLIPSE_T(adj_ind,:);

            % Make profile that includes NEPCs 1-5
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

            for k = 1:5
                if k == 5
                    NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                else
                    NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                end
            end

            R_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            R_int = normintV2(R_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

            % Check shape
            if plotting == 1
                polarplot(Tnew,R_ext,'r');
                hold on
                polarplot(Tnew,R_int,'r');
                polarplot(Tnew,base_ext,'b');
                polarplot(Tnew,base_int,'b');
                title('Changing NEPC 5');
                pause();
                close;
            end

            % Create cases
            case_num = case_num + 1;
            Script = Template; % Reset the script template    
            make_case(case_num,adj_ind,ID,GROUP,base_ext,base_int,Tnew,Script,Erind,Epith);
        
        end
    end

    group = group + 1;
end

set(0,'DefaultFigureWindowStyle','normal');



end




%% Local functions %%
function make_case(case_num,i,ID,GROUP,R_ext,R_int,T,Script,Erind,Epith)
% FILENAME: make_case.m
% AUTHOR: Ryan Larson
% DATE: 11/25/19
%
% PURPOSE: Make a customized Python script corresponding to a unique model.
%       Takes a Python script template defined by a script like
%       write_Python_template3.m or write_Python_template4.m and edits it
%       according to the data for the input model (geometric profile
%       information, material property information, names, etc.). The
%       script template is a column vector of cell arrays that each contain
%       a string, corresponding to a row of the resulting Python script. 
% 
% 
% INPUTS:
%       case_num: Case number (model approximation number)
% 
%       i: Adjusted index of the current cross-section
% 
%       ID: String version of the stalk number
% 
%       GROUP: String version of the group number (corresponds to the index
%       of the current slice location in the slice locations vector)
% 
%       R_ext: Exterior boundary vector
% 
%       R_int: Interior boundary vector
% 
%       T: Theta vector
% 
%       Script: Template script in cell array form
% 
%       Erind: Rind modulus of elasticity
% 
%       Epith: Pith modulus of elasticity
%       
% OUTPUTS:
%        - Creates a customized Python script that can be directly run in
%        Abaqus
%
% NOTES:
%      
% 
% 
% -------------------------------------------------------------------------
% SUBROUTINES:
%   writespline_V2.m: Take an array of XY data, where X is column 1 and Y
%   is column 2, and turn those data points into a list of ordered pairs.
%   This results in one very long string that can be inserted into the
%   Python script; Abaqus uses this to define a spline profile when
%   building a shape. 
% 
% 
% PSEUDO-CODE:
%   Make string versions of the case number, job name, and script name.
% 
%   Convert data from polar to Cartesian coordinates. The resulting Python
%   script will be fed into Abaqus, which expects Cartesian coordinates.
% 
%   Multiply profile values by 1000 to scale to micrometers from
%   millimeters (this is necessary for the transverse models as determined
%   by a mesh convergence study).
% 
%   Transpose XY data to column vectors and combine. X is first column, Y
%   is second column.
% 
%   Repeat the initial data point for exterior and interior to close the
%   profile (all work before this has not repeated the initial point).
% 
%   Determine the data points that should be used for reference points.
%   These should be at the top and bottom of the stalk cross-section when
%   the major diameter of the cross-section is oriented along the X-axis. A
%   simplification is made, where the points closest to theta = 90 degrees
%   and theta = 270 degrees are chosen (this provides flexibility if data
%   sampling is sparse or just doesn't line up exactly on integer degree
%   values).
% 
%   Convert the reference point values to strings for insertion in the
%   Python script.
% 
%   Write the exterior profile as a spline and save as a string.
%   Write the interior profile as a spline and save as a string.
% 
%   Get string versions of rind and pith properties for inserting into the 
%   Python script.
% 
%   Insert all the string variables that have been created into their
%   appropriate positions in the Python script template (still a cell
%   array).
% 
%   Turn the cell array into an actual Python script with a unique name and
%   save. The file will be produced in the current working folder.
% 
% -------------------------------------------------------------------------
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
    % Get string versions of case number, job name, and script name
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

    % Scale units to micrometers from millimeters (allows the mesh size to
    % get small enough based on the mesh convergence study)
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
    
    % Convert the reference point values to strings
    RP1X = sprintf('%0.5g',X_ext(ind90));
    RP1Y = sprintf('%0.5g',Y_ext(ind90));
    RP2X = sprintf('%0.5g',X_ext(ind270));
    RP2Y = sprintf('%0.5g',Y_ext(ind270));

    % Write the spline points and save as a string
    outer_spline = writespline_V2(section_ext);
    inner_spline = writespline_V2(section_int);
    
    % Get string versions of rind and pith properties for inserting into
    % the Python script
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



function [spline] = writespline_V2(data)
% FILENAME: writespline_V2.m
% AUTHOR: Ryan Larson
% DATE: 5/29/19
%
% PURPOSE: Turn spline data from a 2-column array into a string
% 
% 
% INPUTS:
%       data: The original boundary array (2 columns, where column 1 is x
%       data and column 2 is y data)
%       
% OUTPUTS:
%       spline: A string version of data that can be inserted in a Python
%       script
%
% NOTES:
%      
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Create an empty string to hold the spline data (or the ordered XY pairs
%   that define the exterior or interior cross-section profile).
% 
%   For each row in the data (a 2-column array of XY data):
%       Concatenate the current data point to the spline in the ordered
%       pair form, like (x,y).
%   end
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - Writes spline to a text file which can then be copied manually into 
% V2 - Made writespline a function that works with existing functions
% instead of writing the spline to a text file
% V3 - 
%
% -------------------------------------------------------------------------
    %define empty spline and number of x-y points
    spline = '';
    
    %run through 1-column arrays of the x and y data points for the spline, and add to the end of the string with the correct formatting
    for i = 1:length(data)
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), ');
    end
end

function [r] = rpts(N,theta,dmaj,dmin)
% FILENAME: rpts.m
% AUTHOR: Ryan Larson
% DATE: 12/4/19
%
% PURPOSE: Simple function to turn ellipse diameter parameters into polar
% values, given a desired number of points.
% 
% 
% INPUTS:
%       N: The number of points in theta (redundant but haven't gone back
%       to change this in everything downstream)
% 
%       theta: Linearly-spaced vector of theta values from 0 to 2*pi.
% 
%       dmaj: Major diameter of the ellipse.
% 
%       dmin: Minor diameter of the ellipse.
%       
% OUTPUTS:
%       r: The radial value of the ellipse at each value of theta.
%
% NOTES:
%       
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Create an empty row vector of radius values with length N.
% 
%   for each element of the empty row vector:
%       Calculate the radial value of the ellipse, given the major and
%       minor diameters and the current angle.
%   end
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end


function [Erind,Epith] = get_materials(method)
% FILENAME: writespline_V2.m
% AUTHOR: Ryan Larson
% DATE: 5/29/19
%
% PURPOSE: Converts from polar to Cartesian
% 
% 
% INPUTS:
%       method: A string to determine the material selection method.
%       Options include: 
%           'random':   Random material properties
%           'min':      Minimum rind, minimum pith
%           'max':      Maximum rind, maximum pith
%           'minpith':  Random rind, minimum pith
%           'maxpith':  Random rind, maximum pith
%           'minrind':  Minimum rind, random pith
%           'maxrind':  Maximum rind, random pith
%           'avg':      Mean rind, mean pith
% 
% OUTPUTS:
%       Erind: Rind modulus
%       
%       Epith: Pith modulus
%
% NOTES:
%      
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Define mean and standard deviation values for rind stiffness (based on
%   Stubbs 2019 values, in units of N/micrometer^2).
% 
%   Define mean and standard deviation values for pith stiffness (based on
%   Stubbs 2019 values, in units of N/micrometer^2).
% 
%   
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
    % Calculate the random material properties from a normal distribution.
    % Bound with 95% confidence interval, calculated from transverse
    % material properties used in another paper (Stubbs 2019).
    Erind_mean = 8.0747e-04; % THESE VALUES ARE IN N/micrometer^2
    Erind_stdev = 3.3517e-04;
    Erind_95 = [6.7414e-04 9.4081e-04];
    Epith_mean = 2.5976e-05;
    Epith_stdev = 1.0303e-05;
    Epith_95 = [2.1878e-05 3.0075e-05];
    
    % Choose which method to use by the input string
    switch method
        % Generate fully random material properties (both rind and pith
        % random) using the mean and standard deviations for rind and pith.
        % Make sure that the values generated are within the 95% confidence
        % interval.
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
        
        % Use lower bound values for both rind and pith
        case 'min'
            Erind = Erind_95(1);
            Epith = Epith_95(1);
            
        % Use upper bound values for both rind and pith
        case 'max'
            Erind = Erind_95(2);
            Epith = Epith_95(2);
            
        % Use lower bound value for rind and upper bound value for pith
        case 'minrind_maxpith'
            Erind = Erind_95(1);
            Epith = Epith_95(2);
            
        % Use upper bound value for rind and lower bound value for pith
        case 'maxrind_minpith'
            Erind = Erind_95(2);
            Epith = Epith_95(1);
            
        % Use lower bound value for pith and mean value for rind
        case 'minpith'
            Erind = Erind_mean;
            Epith = Epith_95(1);
        
        % Use upper bound value for pith and mean value for rind
        case 'maxpith'
            Erind = Erind_mean;
            Epith = Epith_95(2);
        
        % Use lower bound value for rind and mean value for pith
        case 'minrind'
            Erind = Erind_95(1);
            Epith = Epith_mean;
            
        % Use upper bound value for rind and mean value for pith
        case 'maxrind'
            Erind = Erind_95(2);
            Epith = Epith_mean;
        
        % Use mean value for rind and mean value for pith
        case 'avg'
            Erind = Erind_mean;
            Epith = Epith_mean;
    end
    
end