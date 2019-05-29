function make_case(case_num,i,ID,R_ext,R_int,T,Script)
    CASE = sprintf('%d',case_num);
    jobname = strcat('''Section_',ID,'_',CASE,'''');
    scriptname = strcat('Section_',ID,'_',CASE,'.py');
    
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

    % Transpose data and combine xy
    section_ext = [X_ext', Y_ext'];
    section_int = [X_int', Y_int'];

    % Repeat the last points to close the loop
    section_ext = [section_ext; section_ext(1,:)];
    section_int = [section_int; section_int(1,:)];

    % Get the reference point values in Cartesian coordinates for
    % reference points at 90 and 270 degrees (indices 91 and 271)
    RP1X = sprintf('%0.5g',X_ext(91));
    RP1Y = sprintf('%0.5g',Y_ext(91));
    RP2X = sprintf('%0.5g',X_ext(271));
    RP2Y = sprintf('%0.5g',Y_ext(271));

    % Write the spline points and save as a string
    S = size(section_ext);
    len = S(1);
    outer_spline = writespline_V2(len,section_ext);
    inner_spline = writespline_V2(len,section_int);

    % Insert the case-specific values into the appropriate parts of the
    % Python script template (must be strings)
    Script(17,1) = strcat(Script(17,1),jobname);
    Script(21,1) = strcat(Script(21,1),ID);
    Script(23,1) = strcat(Script(23,1),CASE);
    Script(35,1) = strcat(Script(35,1),RP1X);
    Script(37,1) = strcat(Script(37,1),RP1Y);
    Script(39,1) = strcat(Script(39,1),RP2X);
    Script(41,1) = strcat(Script(41,1),RP2Y);
    Script(61,1) = strcat(Script(61,1),outer_spline);
    Script(84,1) = strcat(Script(84,1),inner_spline);
    
    % Write Python script from the cell array
    filePh = fopen(scriptname,'w');
    fprintf(filePh,'%s\n',Script{:});
    fclose(filePh);
    
end