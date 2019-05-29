function ellipse_sensitivityV2(FileName,numsections,startsection)
% Take in a mat file with ellipse fit data (such as fivefits_top1.mat) and
% create new data points that represent a shift in either a, b, or t, as
% well as the original data points. Follows the organization and methods
% found in create_cases.m.


% numsections: How many cross-sections to create data for (should be same
% value as in create_cases.m)
% startsection: Index of the cross-section to start with (should be same
% value as in create_cases.m)

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'A','B','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int','AVG_RIND_T');

endsection = startsection + numsections - 1;

percent_change = 0.05;
plus_change = 1 + percent_change;
minus_change = 1 - percent_change;

write_Python_template;  % Create Template cell array that can be copied and used to make individualized Python scripts

%% Create all geometry cases for a given cross section
    % Step through the cross sections
    for i = startsection:endsection
        ID = sprintf('%d',i); % Cross-section number

        %% Pure ellipse fit (case 0)
        case_num = 0;
        Script = Template; % Reset the script template
        make_case(case_num,i,ID,ELLIPSE_R_ext,ELLIPSE_R_int,ELLIPSE_T,Script)

        %% Change A (cases 1 and 2: +/-)
        % Adjust exterior points in polar
        Tnew = ELLIPSE_T(i,:);
        Aplus_ext      = rpts(360,ELLIPSE_T(i,:),(plus_change*A(i)),B(i));
        Aminus_ext     = rpts(360,ELLIPSE_T(i,:),(minus_change*A(i)),B(i));

        % Calculate the interior points
        Aplus_int      = Aplus_ext - AVG_RIND_T(i);
        Aminus_int     = Aminus_ext - AVG_RIND_T(i);


%         % Check shape
%         polarplot(Tnew,Aplus_ext);
%         hold on
%         polarplot(Tnew,Aplus_int);
%         pause();
%         polarplot(Tnew,Aminus_ext);
%         polarplot(Tnew,Aminus_int);
%         pause();
%         close;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,Aplus_ext,Aplus_int,Tnew,Script);

        case_num = case_num + 1;
        Script = Template; % Reset the script template
        make_case(case_num,i,ID,Aminus_ext,Aminus_int,Tnew,Script);


        %% Change B (cases 3 and 4: +/-)
        % Adjust exterior points in polar
        Tnew = ELLIPSE_T(i,:);
        Bplus_ext      = rpts(360,ELLIPSE_T(i,:),(A(i)),plus_change*B(i));
        Bminus_ext     = rpts(360,ELLIPSE_T(i,:),(A(i)),minus_change*B(i));

        % Calculate the interior points
        Bplus_int      = Bplus_ext - AVG_RIND_T(i);
        Bminus_int     = Bminus_ext - AVG_RIND_T(i);


%         % Check shape
%         polarplot(Tnew,Bplus_ext);
%         hold on
%         polarplot(Tnew,Bplus_int);
%         pause();
%         polarplot(Tnew,Bminus_ext);
%         polarplot(Tnew,Bminus_int);
%         pause();
%         close;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,Bplus_ext,Bplus_int,Tnew,Script);

        case_num = case_num + 1;
        Script = Template; % Reset the script template
        make_case(case_num,i,ID,Bminus_ext,Bminus_int,Tnew,Script);


        %% Change T (cases 5 and 6: +/-)
        Tnew = ELLIPSE_T(i,:);
        base_ext = rpts(360,ELLIPSE_T(i,:),(A(i)),B(i));

        % Calculate the interior points
        Tplus_int      = base_ext - plus_change*AVG_RIND_T(i);
        Tminus_int     = base_ext - minus_change*AVG_RIND_T(i);

%         % Check shape
%         polarplot(Tnew,base_ext);
%         hold on
%         polarplot(Tnew,Tplus_int);
%         pause();
%         polarplot(Tnew,Tminus_int);
%         pause();
%         close;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,base_ext,Tplus_int,Tnew,Script);

        case_num = case_num + 1;
        Script = Template; % Reset the script template
        make_case(case_num,i,ID,base_ext,Tminus_int,Tnew,Script);


    end


end




%% Local functions %%
function make_case(case_num,i,ID,R_ext,R_int,T,Script)
    CASE = sprintf('%d',case_num);
    jobname = strcat('''Sensitivity_',ID,'_',CASE,'''');
    scriptname = strcat('Sensitivity_',ID,'_',CASE,'.py');
    
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

function [r] = rpts(N,theta,dmaj,dmin)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end
