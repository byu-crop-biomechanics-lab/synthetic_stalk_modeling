function [problem_slice_stalk] = transverse_wrapper_V4(slices,stalknums,material_method,AllSlicesPCA)
% FILENAME: transverse_wrapper_V4.m
% AUTHOR: Ryan Larson
% DATE: 1/17/2020
%
% PURPOSE: Take a subset of the slice locations used in PCA and generate
% all the desired approximations. Creates Python scripts that can be fed
% into Abaqus directly or run as a batch using MasterScript3.py.
% 
% 
% INPUTS:
%       slices: Row vector of slice locations, relative to the node. This 
%       should be a subset of the input that went into AllTransversePCA.m
%       (as of 1/22/2020, this was
%       [-40 -30 -20 -15 -10 -5 0 5 10 15 20 30 40], so slices could be
%       something like [-20 -10 0 30] if used with AllTransversePCA.m)
%       
%       stalknums: A vector of unique integers from 1 to 980 that determines
%       which stalks to sample from (use randperm(980,K) to choose K
%       unique integers from 1 to 980)
%
%       material_method: A string to call out what material selection
%       method to use. Options include:
%           'random' - Randomly sample rind and pith properties
%           'min' - Minimum rind and minimum pith properties
%           'max' - Maximum rind and minimum pith properties
%           'minpith' - Minimum pith properties, random rind properties
%           'maxpith' - Maximum pith properties, random rind properties
%           'minrind' - Minimum rind properties, random pith properties
%           'maxrind' - Maximum rind properties, random pith properties
%           'avg' - Mean rind properties, mean pith properties
% 
%       AllSlicesPCA: AllSlicesPCA.mat
%       
% OUTPUTS:
%       problem_slice_stalk: Indices of stalk slices that had issues during
%       the process. Only useful for debugging.
% 
% NOTES: 
%       Other outputs:
%           - Several .mat files with variables saved from the steps in the
%           process. These are made available for troubleshooting purposes.
%           - Lots of Python scripts, corresponding to the stalks used and 
%           the cases examined. The names do not contain the slice distance
%           info. 11 scripts are created per stalk.
% 
% -------------------------------------------------------------------------
% SUBROUTINES:
%       make_case.m: Create a specimen-specific, model-specific Python
%       script that runs the chosen model in transverse compression when
%       fed to ABAQUS.
% 
%       writespline_V2.m: A subroutine of make_case.m. Convert boundary
%       data from Matlab arrays to strings that can be used in the Python
%       scripts.
% 
%       get_materials.m: Generate rind and pith stiffnesses to use in a
%       given model.
% 
% PSEUDO-CODE:
%   Choose the number of principal components to include in model creation.
%   Load existing PCA data.
% 
%   for each chosen slice distance:
%       Find the index in slice_dists where the current slice exists (this
%       gives indexing information that helps find specific
%       cross-sections).
%       
%       Find the starting index in the PCA data for the current slice location
%       
%       for each stalk:
%           Determine the index of the current slice and stalk in the PCA
%           data.
% 
%           Write string variables for the group number (to tell slice
%           distances apart from each other), the stalk number.
% 
%           Create a copy of the appropriate Python script template.
% 
%           Generate the material properties to use for the current
%           specimen-specific cross-section.
% 
%           Generate the Python script for the real cross-section boundary
%           (case 0).
% 
%           Generate the Python script for the ellipse fit of the current
%           cross-section (case 1).
% 
%           for j = 1:(number of principal components to include)
%               Generate the Python script for the current case, where the
%               ellipse and j principal components are combined.
%           end
% 
%           for j = 2:(number of principal components to include)
%               Generate the Python script for the current case, where the
%               ellipse and the jth principal component are combined.
%           end
%       
%       Change the group number when moving on to the next slice location.
% 
%       end
% 
%   end
%       
% -------------------------------------------------------------------------
% 
% VERSION HISTORY:
% V1 - Choose a continuous range of stalk numbers to sample
% V2 - Input a vector of unique random integers to determine which stalks
% to sample
% V3 - Added group input to number the output scripts for putting in a
% large folder of scripts for different purposes
% V4 - Changed the type of data that gets fed into transverse_wrapper
%
% -------------------------------------------------------------------------

%% Initial variables
% Groups are numbered starting at 1, and are a convenient way to keep data
% separate. Each slice location chosen for analysis initiates a new group.
group = 1;
% Number of principal components to use for analysis.
numNEPCs = 5;
load(AllSlicesPCA);

% Initialize holding variable for slice indices that cause problems
problem_slice_stalk = [];

set(0,'DefaultFigureWindowStyle','docked');

%% Process
% Iterate through slices (determine group number here)
for slice = slices
    
    % Determine the indices in the data where the slice location lives
    sliceidx = find(slice_dists == slice);
    
    % slice_startstop is a variable that loads from AllSlicesPCA.mat.
    % startidx is the index 
    startidx = slice_startstop(sliceidx,2);
    
    % For each slice position, iterate through stalknums
    for stalk = stalknums
        % Get the actual index of the chosen data and create a Python script for
        % that case, numbering by group
        indices = cell2mat(adj_indices(sliceidx,1));
        stalkidx = find(indices == stalk);
        
        if isempty(stalkidx)
            problem_slice_stalk = [problem_slice_stalk; slice, stalk];
            continue
        end
        
        adj_ind = startidx + stalkidx - 1;
        
        %% Create case from ellipse and PCA data (using "ALL" variables)
        
        GROUP = sprintf('%d',group); % Group number
        ID = sprintf('%d',stalk); % Cross-section number
        
        % Create an instance of a blank Python script template in cell
        % array form
        write_Python_template3;
        
        % Select rind and pith properties from a distribution, using the
        % chosen method
        [Erind,Epith] = get_materials(material_method);

        % Real cross section (case 0)
        case_num = 0; % increment this for each case within each cross section
        % Make an instance of Template that will be modified for each
        % unique Python script
        Script = Template;
        % Make the unique Python script for the current case and
        % cross-section
        make_case(case_num,adj_ind,ID,GROUP,ALL_R_ext,ALL_R_int,ALL_ELLIPSE_T,Script,Erind,Epith);

        % Pure ellipse fit (case 1)
        case_num = case_num + 1;
        Script = Template; % Reset the script template
        make_case(case_num,adj_ind,ID,GROUP,ALL_ELLIPSE_R_ext,ALL_ELLIPSE_R_int,ALL_ELLIPSE_T,Script,Erind,Epith);

        % Combined PC cases
        for j = 1:numNEPCs
            case_num = case_num + 1;
            Script = Template; % Reset the script template

            % Calculate the cases with PCs cumulatively added into the
            % ellipse fit
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
            for k = 1:j
                % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
                NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
            end
            
            % Construct the new exterior and interior boundaries
            Rnew_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            Rnew_int = normintV2(Rnew_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));
            
            make_case(case_num,adj_ind,ID,GROUP,Rnew_ext,Rnew_int,ALL_ELLIPSE_T,Script,Erind,Epith);

        end


        % Remaining individual NEPC cases
        for j = 2:numNEPCs
            case_num = case_num + 1;
            Script = Template; % Reset the script template

            % Add the current NEPC to the ellipse in polar coordinates
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
            NEPC_ext = ext_rhocoeffs(adj_ind,j)*ext_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
            
            Rnew_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            Rnew_int = normintV2(Rnew_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));
            make_case(case_num,adj_ind,ID,GROUP,Rnew_ext,Rnew_int,ALL_ELLIPSE_T,Script,Erind,Epith);

        end
        
        
        
        
        
        
    end
    
    % Increase the group number when the slice location changes. This is to
    % prevent multiple scripts with the same name overwriting each other.
    group = group + 1;
end

set(0,'DefaultFigureWindowStyle','normal');

end



%% Localizing all functions used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


function [spline] = writespline_V2(len,data)
% FILENAME: writespline_V2.m
% AUTHOR: Ryan Larson
% DATE: 5/29/19
%
% PURPOSE: Turn vector spline data into a string
% 
% 
% INPUTS:
%       len: Length of the data (a holdover from previous versions, and
%       isn't fully necessary for good code)
% 
%       data: The original boundary vector
%       
% OUTPUTS:
%       spline: A string version of data that can be inserted in a Python
%       script
%
% NOTES:
%      
% 
% 
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
    for i = 1:len 
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), '); 
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