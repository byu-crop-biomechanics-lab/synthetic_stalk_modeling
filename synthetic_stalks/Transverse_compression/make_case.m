function make_case(case_num,i,ID,GROUP,R_ext,R_int,T,Script,Erind,Epith)
% FILENAME: make_case.m
% AUTHOR: Ryan Larson
% DATE: 11/25/19
%
% PURPOSE: Make a customized Python script corresponding to a unique model
% 
% 
% INPUTS:
%       case_num: Case number (model approximation number)
% 
%       i: Adjusted index of the cross-section being worked with
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
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
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
    S = size(section_ext);
    len = S(1);
    outer_spline = writespline_V2(len,section_ext);
    inner_spline = writespline_V2(len,section_int);
    
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