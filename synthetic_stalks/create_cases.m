% create_cases.m: Calculate the necessary information to include in the
% Python script template 

clear; close;
load NEPCs_bottom_987.mat
load Ellipse_fits_bottom_987.mat

numsections = 50;   % Choose the number of cross sections to examine
startsection = 1;   % Choose the starting index of the 50 cross sections
endsection = startsection + numsections - 1;
numNEPCs = 5;       % Choose the number of NEPCs to use in case creation (5 was the original choice)

write_Python_template;  % Create Template cell array that can be copied and used to make individualized Python scripts

%% Create all geometry cases for a given cross section
% Step through the cross sections
for i = startsection:endsection
    ID = sprintf('%d',i); % Cross-section number
    
    %% Real cross section (case 0)
    case_num = 0; % increment this for each case within each cross section
    Script = Template;
    make_case(case_num,i,ID,R_ext,R_int,ELLIPSE_T,Script)
    
    %% Pure ellipse fit (case 1)
    case_num = case_num + 1;
    Script = Template; % Reset the script template
    jobname = strcat('Section_',ID,'_',case_num);
    make_case(case_num,i,ID,ELLIPSE_R_ext,ELLIPSE_R_int,ELLIPSE_T,Script)
    
    
    %% Combined NEPC cases
    for j = 1:numNEPCs
        case_num = case_num + 1;
        Script = Template; % Reset the script template
        
        % Calculate the cases with NEPCs cumulatively added into the
        % ellipse fit
        NEPC = zeros(1,size(ext_rhoPCAs,1));
        for k = 1:j
            % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
            NEPC = NEPC + ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
            
        end
        
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);

        make_case(case_num,i,ID,Rnew_ext,Rnew_int,ELLIPSE_T,Script)
        
    end
    
    
    %% Remaining individual NEPC cases
    for j = 2:numNEPCs
        case_num = case_num + 1;
        Script = Template; % Reset the script template
        
        % Add the current NEPC to the ellipse in polar coordinates
        NEPC = zeros(1,size(ext_rhoPCAs,1));
        NEPC = ext_rhocoeffs(i,j)*ext_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        make_case(case_num,i,ID,Rnew_ext,Rnew_int,ELLIPSE_T,Script)
        
    end
    
    
    
    
end


%% Local functions %%
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


% 
% %% Individual effects of NEPCs
% % Step through the cross sections
% for i = 1:length(indices)
%     NEPC = zeros(1,size(ext_rhoPCAs,1));
%     % Step through the NEPCs
%     for j = 1:size(ext_rhocoeffs,1)
%         
%         % Add the current NEPC to the ellipse in polar coordinates
%         NEPC = ext_rhocoeffs(i,j)*ext_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
%         Rnew_ext = ELLIPSE_R(i,:) - NEPC;
%         Rnew_int = Rnew_ext - AVG_RIND_T(i);
%         
%         % Convert from polar to Cartesian
%         xy_ext = convert_to_xy(Rnew_ext',ELLIPSE_T(i,:)');
%         xy_int = convert_to_xy(Rnew_int',ELLIPSE_T(i,:)');
%         xy_original = [ELLIPSE_XY(i,:,1)',ELLIPSE_XY(i,:,2)'];
% 
%         % Repeat the last points to close the loop
%         xy_ext   = [xy_ext; xy_ext(1,:)];
%         xy_int   = [xy_int; xy_int(1,:)];
%         xy_original = [xy_original; xy_original(1,:)];
%         
%         % Plot to check the result
%         plot(xy_ext(:,1),xy_ext(:,2));
%         hold on
%         plot(xy_int(:,1),xy_int(:,2));
%         plot(xy_original(:,1),xy_original(:,2),'m');
%         axis equal
%         titlestring = sprintf('Ellipse %d with PC%d',i,j);
%         title(titlestring);
%         pause();
%         hold off
%         
% %         % Export a .txt file with the new shape data for exterior and
% %         % interior boundary points
% %         % Save base points as txt file
% %         S = size(xy_ext);
% %         len = max(S);
% %         outputName_ext = sprintf('Ellipse_wPC%d_ext',j);
% %         outputName_int = sprintf('Ellipse_wPC%d_int',j);
% %         % (also export a mat file with the same data for finding the
% %         % reference points to use in Abaqus)
% %         matname_ext = strcat(outputName_ext,num2str(i),'.mat');
% %         matname_int = strcat(outputName_int,num2str(i),'.mat');
% %         save(matname_ext,'xy_ext');
% %         save(matname_int,'xy_int');
% %         writespline(len,xy_ext,i,outputName_ext)
% %         writespline(len,xy_int,i,outputName_int)
%     end
%     
%     
% end
% 
% 
% %% Convergence effects of NEPCs
% 
% % Step through the cross sections
% for i = 1:length(indices)
%     NEPC = zeros(1,size(ext_rhoPCAs,1));
%     % Step through the NEPCs
%     for j = 1:size(ext_rhocoeffs,1)
%         
%         % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
%         NEPC = NEPC + ext_rhocoeffs(i,j)*ext_rhoPCAs(:,j)';
%         Rnew_ext = ELLIPSE_R(i,:) - NEPC;
%         Rnew_int = Rnew_ext - AVG_RIND_T(i);
%         
%         % Convert from polar to Cartesian
%         xy_ext = convert_to_xy(Rnew_ext',ELLIPSE_T(i,:)');
%         xy_int = convert_to_xy(Rnew_int',ELLIPSE_T(i,:)');
%         xy_original = [ELLIPSE_XY(i,:,1)',ELLIPSE_XY(i,:,2)'];
% 
%         % Repeat the last points to close the loop
%         xy_ext   = [xy_ext; xy_ext(1,:)];
%         xy_int   = [xy_int; xy_int(1,:)];
%         xy_original = [xy_original; xy_original(1,:)];
%         
%         % Plot to check the result
%         plot(xy_ext(:,1),xy_ext(:,2));
%         hold on
%         plot(xy_int(:,1),xy_int(:,2));
%         plot(xy_original(:,1),xy_original(:,2),'m');
%         axis equal
%         titlestring = sprintf('Ellipse %d with all PCs through PC%d',i,j);
%         title(titlestring);
%         pause();
%         hold off
%         
% %         % Export a .txt file with the new shape data for exterior and
% %         % interior boundary points
% %         % Save base points as txt file
% %         S = size(xy_ext);
% %         len = max(S);
% %         outputName_ext = sprintf('Ellipse_wPC1_to_%d_ext',j);
% %         outputName_int = sprintf('Ellipse_wPC1_to_%d_int',j);
% %         % (also export a mat file with the same data for finding the
% %         % reference points to use in Abaqus)
% %         matname_ext = strcat(outputName_ext,num2str(i),'.mat');
% %         matname_int = strcat(outputName_int,num2str(i),'.mat');
% %         save(matname_ext,'xy_ext');
% %         save(matname_int,'xy_int');
% %         writespline(len,xy_ext,i,outputName_ext)
% %         writespline(len,xy_int,i,outputName_int)
%     end
%     
%     
% end
% 
% close;