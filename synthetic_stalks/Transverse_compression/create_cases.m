function create_cases(NEPCdata,GoodEllipseData,SelectedData,problem_indices,numNEPCs,material_method,SaveName)
    % create_cases.m: Calculate the necessary information to include in the
    % Python scripts
    
    load(NEPCdata);
    load(GoodEllipseData);
    load(SelectedData);

    N = size(ELLIPSE_T,1);
    MaterialProps = zeros(N,(2 + 2*numNEPCs - 1),2);
    

    write_Python_template;  % Create Template cell array that can be copied and used to make individualized Python scripts

    %% Create all geometry cases for a given cross section
    % Step through the cross sections
    stalks = selectedTable.StkNum;
    
    % Remove the stalk numbers that had ellipse fit problems
    stalks(problem_indices) = [];
    
    for i = 1:N
        
        ID = sprintf('%d',stalks(i)); % Cross-section number

        %% Real cross section (case 0)
        case_num = 0; % increment this for each case within each cross section
        Script = Template;
        [Erind,Epith] = get_materials(material_method);
        make_case(case_num,i,ID,R_ext,R_int,ELLIPSE_T,Script,Erind,Epith);
        MaterialProps(i,case_num+1,1) = Erind;
        MaterialProps(i,case_num+1,2) = Epith;

        %% Pure ellipse fit (case 1)
        case_num = case_num + 1;
        Script = Template; % Reset the script template
        [Erind,Epith] = get_materials(material_method);
        make_case(case_num,i,ID,ELLIPSE_R_ext,ELLIPSE_R_int,ELLIPSE_T,Script,Erind,Epith);
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
                NEPC_ext = NEPC_ext + ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
%                 NEPC_int = NEPC_int + int_rhocoeffs(i,k)*int_rhoPCAs(:,k)';
            end

            Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
            Rnew_int = Rnew_ext - AVG_RIND_T(i);
            
            [Erind,Epith] = get_materials(material_method);
            make_case(case_num,i,ID,Rnew_ext,Rnew_int,ELLIPSE_T,Script,Erind,Epith);
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
            NEPC_ext = ext_rhocoeffs(i,j)*ext_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
%             NEPC_int = int_rhocoeffs(i,j)*int_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
            Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
%             Rnew_int = Rnew_ext - AVG_RIND_T(i);

            [Erind,Epith] = get_materials(material_method);
            make_case(case_num,i,ID,Rnew_ext,Rnew_int,ELLIPSE_T,Script,Erind,Epith);
            MaterialProps(i,case_num+1,1) = Erind;
            MaterialProps(i,case_num+1,2) = Epith;

        end

    end
    
    % Save the final data in a new mat file
    FolderName = pwd;
    SaveFile = fullfile(FolderName, SaveName);
    save(SaveFile,'MaterialProps');
    
    

end
    
%% Local functions %%
function make_case(case_num,i,ID,R_ext,R_int,T,Script,Erind,Epith)
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
    Script(21,1) = strcat(Script(21,1),ID);
    Script(23,1) = strcat(Script(23,1),CASE);
    Script(31,1) = strcat(Script(31,1),rindE);
    Script(33,1) = strcat(Script(33,1),pithE);
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
    
    end
    
end