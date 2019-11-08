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