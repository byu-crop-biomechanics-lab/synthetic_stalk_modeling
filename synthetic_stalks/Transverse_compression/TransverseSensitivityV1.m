function TransverseSensitivityV1(slices,stalknums,AllSlicesPCA,percent_change)
% FILENAME: TransverseSensitivityV1.m
% AUTHOR: Ryan Larson
% DATE: 1/24/2020
%
% PURPOSE: Wrap the majority of the data production process into a single
% script
% 
% 
% INPUTS:
%       slices - A subset of the input that went into AllTransversePCA.m (as of
%       1/22/2020, this was [-40 -30 -20 -15 -10 -5 0 5 10 15 20 30 40])
%       stalknums - A vector of unique integers from 1 to 980 that determines
%       which stalks to sample from (use randperm(980,K) to choose K
%       unique integers from 1 to 980)
%       AllSlicesPCA - AllSlicesPCA.mat
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
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

%% Initial variables

load(AllSlicesPCA);
group = 1;
numNEPCs = 5;


% percent_change = 0.05;
plus_change = 1 + percent_change;
minus_change = 1 - percent_change;

write_Python_template3;  % Create Template cell array that can be copied and used to make individualized Python scripts

Rind = zeros(length(A),11);
Pith = zeros(size(Rind));

%% Create all geometry cases for a given cross section
    % Step through the cross sections
    for i = 1:length(A)
        GROUP = sprintf('%d',group); % Group number
        ID = sprintf('%d',stalknums(i)); % Cross-section number
        
        % Get material properties
        [Erind,Epith] = get_materials('random');
        
        % Make profile that includes NEPCs 1-5
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:5
            % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
            NEPC_ext = NEPC_ext + ext_rhocoeffs(newgoodstalknums(i),k)*ext_rhoPCAs(:,k)';
        end
        
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);

        
        %% Base case (case 0)
        case_num = 0;
        Script = Template; % Reset the script template
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,ELLIPSE_T,Script,Erind,Epith)
        Rind(i,1) = Erind;
        Pith(i,1) = Epith;

        
        %% Change A (case 1)
        % Adjust exterior points in polar
        Tnew = ELLIPSE_T(i,:);
        Aplus_ellipse = rpts(360,ELLIPSE_T(i,:),(plus_change*A(i)),B(i));
        Aplus_ext = Aplus_ellipse - NEPC_ext;

        % Calculate the interior points
        Aplus_int = Aplus_ext - AVG_RIND_T(i);

        % Check shape
        if plotting == 1
            polarplot(Tnew,Aplus_ext);
            hold on
            polarplot(Tnew,Aplus_int);
            polarplot(Tnew,Rnew_ext);
            polarplot(Tnew,Rnew_int);
            pause();
            close;
        end

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Aplus_ext,Aplus_int,Tnew,Script,Erind,Epith);
        Rind(i,2) = Erind;
        Pith(i,2) = Epith;
        

        %% Change B (case 2)
        % Adjust exterior points in polar
        Tnew = ELLIPSE_T(i,:);
        Bplus_ellipse = rpts(360,ELLIPSE_T(i,:),(A(i)),plus_change*B(i));
        Bplus_ext = Bplus_ellipse - NEPC_ext;

        % Calculate the interior points
        Bplus_int      = Bplus_ext - AVG_RIND_T(i);

        % Check shape
        if plotting == 1
            polarplot(Tnew,Bplus_ext);
            hold on
            polarplot(Tnew,Bplus_int);
            polarplot(Tnew,Rnew_ext);
            polarplot(Tnew,Rnew_int);
            pause();
            close;
        end

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Bplus_ext,Bplus_int,Tnew,Script,Erind,Epith);
        Rind(i,3) = Erind;
        Pith(i,3) = Epith;

        
        %% Change T (case 3)
        Tnew = ELLIPSE_T(i,:);
        base_ext = Rnew_ext;

        % Calculate the interior points
        Tplus_int      = base_ext - plus_change*AVG_RIND_T(i);

        % Check shape
        if plotting == 1
            polarplot(Tnew,base_ext);
            hold on
            polarplot(Tnew,Tplus_int);
            polarplot(Tnew,Rnew_ext);
            polarplot(Tnew,Rnew_int);
            pause();
            close;
        end

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,base_ext,Tplus_int,Tnew,Script,Erind,Epith);
        Rind(i,4) = Erind;
        Pith(i,4) = Epith;
        
        
        %% Change Erind (case 4)
        Tnew = ELLIPSE_T(i,:);
        
        % Calculate the new Erind
        Erind_plus = Erind*plus_change;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,Tnew,Script,Erind_plus,Epith);
        Rind(i,5) = Erind_plus;
        Pith(i,5) = Epith;
        
        
        %% Change Epith (case 5)
        Tnew = ELLIPSE_T(i,:);
        
        % Calculate the new Erind
        Epith_plus = Epith*plus_change;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,Tnew,Script,Erind,Epith_plus);
        Rind(i,6) = Erind;
        Pith(i,6) = Epith_plus;
        
        
        %% Change NEPC 1 (case 6)
        Tnew = ELLIPSE_T(i,:);
        
        % Make profile that includes NEPCs 1-5
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:5
            if k == 1
                NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(newgoodstalknums(i),k)*ext_rhoPCAs(:,k)';
            else
                NEPC_ext = NEPC_ext + ext_rhocoeffs(newgoodstalknums(i),k)*ext_rhoPCAs(:,k)';
            end
        end
        
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,Tnew,Script,Erind,Epith);
        Rind(i,7) = Erind;
        Pith(i,7) = Epith;
        
        
        %% Change NEPC 2 (case 7)
        Tnew = ELLIPSE_T(i,:);
        
        % Make profile that includes NEPCs 1-5
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:5
            if k == 2
                NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(newgoodstalknums(i),k)*ext_rhoPCAs(:,k)';
            else
                NEPC_ext = NEPC_ext + ext_rhocoeffs(newgoodstalknums(i),k)*ext_rhoPCAs(:,k)';
            end
        end
        
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,Tnew,Script,Erind,Epith);
        Rind(i,8) = Erind;
        Pith(i,8) = Epith;
        
        
        %% Change NEPC 3 (case 8)
        Tnew = ELLIPSE_T(i,:);
        
        % Make profile that includes NEPCs 1-5
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:5
            if k == 3
                NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
            else
                NEPC_ext = NEPC_ext + ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
            end
        end
        
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,Tnew,Script,Erind,Epith);
        Rind(i,9) = Erind;
        Pith(i,9) = Epith;
        
        
        %% Change NEPC 4 (case 9)
        Tnew = ELLIPSE_T(i,:);
        
        % Make profile that includes NEPCs 1-5
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:5
            if k == 4
                NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
            else
                NEPC_ext = NEPC_ext + ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
            end
        end
        
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,Tnew,Script,Erind,Epith);
        Rind(i,10) = Erind;
        Pith(i,10) = Epith;
        
        
        %% Change NEPC 5 (case 10)
        Tnew = ELLIPSE_T(i,:);
        
        % Make profile that includes NEPCs 1-5
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:5
            if k == 5
                NEPC_ext = NEPC_ext + plus_change*ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
            else
                NEPC_ext = NEPC_ext + ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
            end
        end
        
        Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,i,ID,GROUP,Rnew_ext,Rnew_int,Tnew,Script,Erind,Epith);
        Rind(i,11) = Erind;
        Pith(i,11) = Epith;
        
        
    end

    FolderName = pwd;
    Materials = 'MaterialsUsed.mat';
    SaveFile = fullfile(FolderName, Materials);
    save(SaveFile,'Rind','Pith');

end




%% Local functions %%
function make_case(case_num,i,ID,GROUP,R_ext,R_int,T,Script,Erind,Epith)
    CASE = sprintf('%d',case_num);
    jobname = strcat('''Group_',GROUP,'_','Sensitivity_',ID,'_',CASE,'''');
    scriptname = strcat('Group_',GROUP,'_','Sensitivity_',ID,'_',CASE,'.py');
    
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
            
        case 'minpith'
            Erind = Erind_mean;
            Epith = Epith_95(1);
            
        case 'maxpith'
            Erind = Erind_mean;
            Epith = Epith_95(2);
            
        case 'minrind'
            Erind = Erind_95(1);
            Epith = Epith_mean;
            
        case 'maxrind'
            Erind = Erind_95(2);
            Epith = Epith_mean;
    
        case 'avg'
            Erind = Erind_mean;
            Epith = Epith_mean;
    end
    
end