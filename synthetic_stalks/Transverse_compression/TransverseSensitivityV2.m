function TransverseSensitivityV2(slices,stalknums,nunique,AllSlicesPCA,percent_change,numNEPCs,plotting)
% FILENAME: TransverseSensitivityV2.m
% AUTHOR: Ryan Larson
% DATE: 3/7/2020
%
% PURPOSE: Updated material sensitivity study
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
%       nunique: The number of unique slices to use in the sensitivity
%       study (try 10 to start - should result in 240 Python scripts)
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
%       rpts.m: Create an ellipse from major and minor diameter values.
% 
%       get_materials.m: Generate rind and pith stiffnesses to use in a
%       given model.
% 
% PSEUDO-CODE:
%   Load PCA data.
% 
%   Randomly sample cross-sections from the data at the slice locations in
%   slices.
%   
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

%% Initial variables
set(0,'DefaultFigureWindowStyle','docked');
load(AllSlicesPCA);

problem_slice_stalk = [];

plus_change = 1 + percent_change;

% slices_unique saves the randomly chosen specific cross-sections (chosen
% from the Group II subset used in FEA)
% Each row is a unique cross-section to be approximated.
slices_unique = zeros(nunique,2);

% The first column is random slice locations, chosen from slices (not from
% the full PCA data, which most likely includes more slice locations than
% the "slices" input to this function)
slices_unique(:,1) = datasample(slices,nunique)';

% The second column is random stalk numbers, chosen from stalknums.
slices_unique(:,2) = datasample(stalknums,nunique)';

% Give numbers to the slice locations used, or groups. These will most
% likely not be 
% groups = [1 2 3 4];
groups = linspace(1,length(slices),length(slices));

write_Python_template4;  % Create Template cell array that can be copied and used to make individualized Python scripts


%% Create all geometry cases for a given cross section
% Iterate through slices (determine group number here)
for i = 1:nunique
    
    % Since the slices for the sensitivity study are randomly generated, we
    % have to access it from slices_unique.
    slice = slices_unique(i,1);
    
    % Convenience variables for finding the desired slice and stalk in the
    % PCA data
    sliceidx = find(slice_dists == slice);
    startidx = slice_startstop(sliceidx,2);
    
    stalk = slices_unique(i,2);
    
    % Get the actual index of the chosen data and create a Python script for
    % that case, numbering by group
    indices = cell2mat(adj_indices(sliceidx,1));
    stalkidx = find(indices == stalk);
    
    % If the randomly-generated cross-section isn't one of the successful
    % ellipse fit cases, then skip it and add the index to
    % problem_slice_stalk so it's clear when that happens.
    if isempty(stalkidx)
        problem_slice_stalk = [problem_slice_stalk; slice, stalk];
        continue
    end
    
    % Calculate the index of the specimen-specific cross-section in the PCA
    % data.
    adj_ind = startidx + stalkidx - 1;
    
    for group = groups
        % Create case from ellipse and PCA data (using "ALL" variables)
        
        GROUP = sprintf('%d',group); % Group number
        ID = sprintf('%d',stalk); % Cross-section number        

        % This is slightly different from write_Python_template3 because
        % the line numbers end up being different for the sensitivity study
        write_Python_template4;

        % Get material properties
        % THIS BLOCK ASSUMES THE FUNCTION IS BEING
        % USED FOR THE MATERIAL SENSITIVITY STUDY, AND NOT A GENERAL
        % SENSITIVITY STUDY FOR ALL PARAMETERS
        if group == 1
            material_method = 'min';
        elseif group == 2
            material_method = 'minrind_maxpith';
        elseif group == 3
            material_method = 'max';
        else
            material_method = 'maxrind_minpith';
        end
        
        [Erind,Epith] = get_materials(material_method);

        % Make profile that includes NEPCs 1-5
        NEPC_ext = zeros(1,size(ext_rhoPCAs,1));

        for k = 1:numNEPCs
            % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
            NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
        end

        base_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
        base_int = normintV2(base_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

    %         Rnew_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
    % %       Rnew_int = ALL_ELLIPSE_R_int(adj_ind,:) - NEPC_ext;
    %         Rnew_int = normintV2(Rnew_ext,ALL_ELLIPSE_T(adj_ind,:),ALL_AVG_RIND_T(adj_ind));

        %% Base case (case 0)
        case_num = 0;
        Script = Template; % Reset the script template
        make_case(case_num,adj_ind,ID,slice,GROUP,base_ext,base_int,ALL_ELLIPSE_T,Script,Erind,Epith)

        %% Change A (case 1)
        % Adjust exterior points in polar
        Tnew = ALL_ELLIPSE_T(adj_ind,:);
        Aplus_ellipse = rpts(360,ALL_ELLIPSE_T(adj_ind,:),(plus_change*ALL_A(adj_ind)),ALL_B(adj_ind));
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
        make_case(case_num,adj_ind,ID,slice,GROUP,Aplus_ext,Aplus_int,Tnew,Script,Erind,Epith);


        %% Change B (case 2)
        % Adjust exterior points in polar
        Tnew = ALL_ELLIPSE_T(adj_ind,:);
        Bplus_ellipse = rpts(360,ALL_ELLIPSE_T(adj_ind,:),(ALL_A(adj_ind)),plus_change*ALL_B(adj_ind));
        Bplus_ext = Bplus_ellipse - NEPC_ext;

        % Calculate the interior points
    %         Bplus_int      = Bplus_ext - AVG_RIND_T(i);
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
        make_case(case_num,adj_ind,ID,slice,GROUP,Bplus_ext,Bplus_int,Tnew,Script,Erind,Epith);


        %% Change T (case 3)
        Tnew = ALL_ELLIPSE_T(adj_ind,:);
    %         base_ext = base_ext;

        % Calculate the interior points
    %         Tplus_int = base_ext - plus_change*ALL_AVG_RIND_T(adj_ind);
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
        make_case(case_num,adj_ind,ID,slice,GROUP,base_ext,Tplus_int,Tnew,Script,Erind,Epith);        

        %% Change Erind (case 4)
        Tnew = ALL_ELLIPSE_T(adj_ind,:);

        % Calculate the new Erind
        Erind_plus = Erind*plus_change;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,adj_ind,ID,slice,GROUP,base_ext,base_int,Tnew,Script,Erind_plus,Epith);        

        %% Change Epith (case 5)
        Tnew = ALL_ELLIPSE_T(adj_ind,:);

        % Calculate the new Erind
        Epith_plus = Epith*plus_change;

        % Create cases
        case_num = case_num + 1;
        Script = Template; % Reset the script template    
        make_case(case_num,adj_ind,ID,slice,GROUP,base_ext,base_int,Tnew,Script,Erind,Epith_plus);        
        
        
        % ADD ADDITIONAL CASES HERE FOR PRINCIPAL COMPONENT
        % PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    end

end



set(0,'DefaultFigureWindowStyle','normal');

FolderName = pwd;
Slices = 'SlicesUsed.mat';
SaveFile = fullfile(FolderName, Slices);
save(SaveFile,'slices_unique');






end




%% Local functions %%
function make_case(case_num,i,ID,slice,GROUP,R_ext,R_int,T,Script,Erind,Epith)
    slicenum = sprintf('%d',slice); % Slice number
    if slice < 0
        slicenum = slicenum(2:end);
        slicenum = strcat('___',slicenum);
        slicenum = sprintf('%s',slicenum);
    end
    CASE = sprintf('%d',case_num);
    jobname = strcat('''Sensitivity_',slicenum,'_',ID,'_','Group_',GROUP,'_',CASE,'''');
    scriptname = strcat('Sensitivity_',slicenum,'_',ID,'_','Group_',GROUP,'_',CASE,'.py');
    
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
%     slicenum = sprintf('''%s''',slicenum);
    slice = sprintf('%d',slice);
    Script(21,1) = strcat(Script(21,1),slice);
    Script(23,1) = strcat(Script(23,1),ID);
    Script(25,1) = strcat(Script(25,1),GROUP);
    Script(27,1) = strcat(Script(27,1),CASE);
    Script(35,1) = strcat(Script(35,1),rindE);
    Script(37,1) = strcat(Script(37,1),pithE);
    Script(39,1) = strcat(Script(39,1),RP1X);
    Script(41,1) = strcat(Script(41,1),RP1Y);
    Script(43,1) = strcat(Script(43,1),RP2X);
    Script(45,1) = strcat(Script(45,1),RP2Y);
    Script(65,1) = strcat(Script(65,1),outer_spline);
    Script(88,1) = strcat(Script(88,1),inner_spline);
    
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
            
        case 'minrind_maxpith'
            Erind = Erind_95(1);
            Epith = Epith_95(2);
            
        case 'maxrind_minpith'
            Erind = Erind_95(2);
            Epith = Epith_95(1);
            
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