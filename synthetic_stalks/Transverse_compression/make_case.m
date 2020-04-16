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


%% Localizing subfunctions
function [spline] = writespline_V2(data)
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
% -------------------------------------------------------------------------
% SUBROUTINES:
%       N/A
% 
% PSEUDO-CODE:
%   Read in data in XY column form (X is column 1, Y is column 2).
% 
%   Create an empty string for the spline data.
% 
%   for each row of data:
%       Add an ordered pair (in string form) to the spline.
%   end
% 
% 
% VERSION HISTORY:
% V1 - Writes spline to a text file which can then be copied manually into
% a Python script (very old original process)
% V2 - Made writespline a function that works with existing functions
% instead of writing the spline to a text file
% V3 - 
%
% -------------------------------------------------------------------------
    %define empty spline and number of x-y points
    spline = '';

    %run through 1-column arrays of the x and y data points for the spline, and add to the end of the string with the correct formatting
    for i = 1:size(data,1)
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), '); 
    end
end