function [] = Bending_Stress(slices,stalknums,numNEPCs,E_ratio)
% FILENAME: Bending_Stress.m
% AUTHOR: Michael Ottesen
% DATE: 5/2020
%
% PURPOSE: Find bending stress percent error for ellipse fit and subsequent
%          principal components.
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
%       numNEPCs: An integer of the desired number of principal components
%       to include in the analysis.
%
%       E_ratio: A number that defines how much greater the rind modulus is
%       compared to the pith modulus. For example, if a modulus ratio Er:Ep
%       of 20:1 is desired, enter '20'.
%       
% OUTPUTS:       
%       S_err_prctile: The values of the key percentiles for the stiffness
%       errors. The percentiles are also listed as [5,25,50,75,90]th
%       percentiles. Each row represents the ellipse with principal
%       components stiffness error as compared to the real shape.
%
% NOTES: 
%       - 
% 
% PSEUDO-CODE:
%   Load AllSicesPCA data.
%   Initialize variables and arrays.
%   Begin large loop.
%       Get cross section data.
%       Find inertia data for ellipse and ellipse plus principal
%       components.
%       Also calculate inertial center in y-direction (for max stress).
%       Determine max distance from center to edge.
%   Calculate stress for all cross sections.
%   Apply modulus ratio and find percent error for data.
%   Create box plot to represent percent error in moment area of inertia of
%   PC cross sections compared to original.
%   Create box plot to represent percent error in max stress of PC cross
%   sections compared to the original.
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------



% Load data from AllTransversePCA.mat
load('AllSlicesPCA.mat')

% Intitialize the problem stalk indices
problem_slice_stalk = [];

% Define slice number to keep trak of number of iterations (slice dists as
% an iteration number)
slicenum = 0;
CSnum = 0;

for slice = slices
    slicenum = slicenum + 1
    
    % Determine the indices in the data where the slice location lives
    sliceidx = find(slice_dists == slice);
    
    % slice_startstop is a variable that loads from AllSlicesPCA.mat.
    % startidx is the index 
    startidx = slice_startstop(sliceidx,2);
    
    % For each slice position, iterate through stalknums
    for stalk = stalknums
        CSnum = CSnum + 1;
        % Get the actual index of the chosen data and create a Python script for
        % that case, numbering by group
        indices = cell2mat(adj_indices(sliceidx,1));
        stalkidx = find(indices == stalk);
        
        if isempty(stalkidx)
            problem_slice_stalk = [problem_slice_stalk; slice, stalk];
            continue
        end
        
        % adj_ind is the row index where the specific cross-section is
        % within the "ALL" arrays
        adj_ind = startidx + stalkidx - 1;
        
        % Real cross section (case 0)
        case_num = 0; % increment this for each case within each cross section
        
        % Track the indice number in the first column
        Ip(CSnum,1,slicenum) = adj_ind;
        Ir(CSnum,1,slicenum) = adj_ind;
        Y(CSnum,1,slicenum) = adj_ind;
        
        
        % Calculate pith and rind inertias for true cross section
        [int_x,int_y] = pol2cart(ALL_ELLIPSE_T(adj_ind,:),ALL_R_int(adj_ind,:));
        [ext_x,ext_y] = pol2cart(ALL_ELLIPSE_T(adj_ind,:),ALL_R_ext(adj_ind,:));
        [~,I_int,~] = polygeom(int_x,int_y);
        [~,I_ext,~] = polygeom(ext_x,ext_y);
        Ip(CSnum,2,slicenum) = I_int(4);
        Ir(CSnum,2,slicenum) = I_ext(4)-I_int(4);
        % Calculate neutral axis for cross section
        [y,~] = Center_Inertia(int_x,int_y,ext_x,ext_y,1e-4);
        y1 = max(ext_y) - y;
        y2 = min(ext_y) - y;
        if abs(y1) > abs(y2)
            Y(CSnum,2,slicenum) = abs(y1);
        elseif abs(y2) > abs(y1)
            Y(CSnum,2,slicenum) = abs(y2);
        end
        % Calculate max bending stress
        
        
        
        % Calculate pith and rind inertias for ellipse cross section
        [int_x,int_y] = pol2cart(ALL_ELLIPSE_T(adj_ind,:),ALL_ELLIPSE_R_int(adj_ind,:));
        [ext_x,ext_y] = pol2cart(ALL_ELLIPSE_T(adj_ind,:),ALL_ELLIPSE_R_ext(adj_ind,:));
        [~,I_int,~] = polygeom(int_x,int_y);
        [~,I_ext,~] = polygeom(ext_x,ext_y);
        Ip(CSnum,3,slicenum) = I_int(4);
        Ir(CSnum,3,slicenum) = I_ext(4)-I_int(4);
        % Calculate neutral axis for cross section
        [y,~] = Center_Inertia(int_x,int_y,ext_x,ext_y,1e-4);
        y1 = max(ext_y) - y;
        y2 = min(ext_y) - y;
        if abs(y1) > abs(y2)
            Y(CSnum,3,slicenum) = abs(y1);
        elseif abs(y2) > abs(y1)
            Y(CSnum,3,slicenum) = abs(y2);
        end

        % Combined PC cases
        for j = 1:numNEPCs
            case_num = case_num + 1;

            % Calculate the cases with PCs cumulatively added into the
            % ellipse fit
            NEPC_int = zeros(1,size(ext_rhoPCAs,1));
            NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
            
            for k = 1:j
                % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
                NEPC_ext = NEPC_ext + ext_rhocoeffs(adj_ind,k)*ext_rhoPCAs(:,k)';
                NEPC_int = NEPC_int + int_rhocoeffs(adj_ind,k)*int_rhoPCAs(:,k)';
            end
            
            % Construct the new exterior and interior boundaries
            Rnew_ext = ALL_ELLIPSE_R_ext(adj_ind,:) - NEPC_ext;
            Rnew_int = ALL_ELLIPSE_R_int(adj_ind,:) - NEPC_int;
            
            % Calculate rind and pith area for adjusted ellipse shape
            [int_x,int_y] = pol2cart(ALL_ELLIPSE_T(adj_ind,:),Rnew_int);
            [ext_x,ext_y] = pol2cart(ALL_ELLIPSE_T(adj_ind,:),Rnew_ext);
            [~,I_int,~] = polygeom(int_x,int_y);
            [~,I_ext,~] = polygeom(ext_x,ext_y);
            Ip(CSnum,j+3,slicenum) = I_int(4);
            Ir(CSnum,j+3,slicenum) = I_ext(4)-I_int(4);
            % Calculate neutral axis for cross section
            [y,~] = Center_Inertia(int_x,int_y,ext_x,ext_y,1e-4);
            y1 = max(ext_y) - y;
            y2 = min(ext_y) - y;
            if abs(y1) > abs(y2)
                Y(CSnum,j+3,slicenum) = abs(y1);
            elseif abs(y2) > abs(y1)
                Y(CSnum,j+3,slicenum) = abs(y2);
            end

            
        end
    end
end

sigma = Y./Ir;

% Apply Modulus ratio to moment area of inertia
ref = E_ratio*Ir(:,2,:) + 1*Ip(:,2,:);
VE = E_ratio*Ir(:,3:numNEPCs+3,:) + 1*Ip(:,3:numNEPCs+3,:);

% Get percent error for area moment of inertia for each varried ellipse 
% compared to real shape
for i = 1:numNEPCs+1
    S_err(:,i,:) =  ((VE(:,i,:) - ref(:,1,:))./ref(:,1,:))*100;
    S_err_prctile_i(i,:,:) = prctile(S_err(:,i,:),[5,25,50,75,95]);
    for j = 1:5
       S_err_prctile(i,j) = mean(S_err_prctile_i(i,j,:)); 
    end
end


% Create labels according to the number of principal components used in
% the study (cumulative cases followed by remaining individual cases)
all_labels = strings(1,(1+numNEPCs));
all_labels(1,1:2) = ["Ellipse","Ellipse + PC 1"];
for i = 2:numNEPCs
    addlabel = "Ellipse + PCs 1-" + num2str(i);
    all_labels(1,i+1) = addlabel;
end

% Combine all percent error data into one matrix for box plotting
S_err_plot = S_err(:,:,1);
for i = 1:length(slices)-1
S_err_plot = [S_err_plot;S_err(:,:,i+1)];
end

% Grab upper and lower limits by using 5th and 95th percentile data
uplimrow = S_err_prctile(:,5)';
lolimrow = S_err_prctile(:,1)';

% Add a buffer between the calculated outer reach of the whiskers and the
% edge of the plot. Round to the nearest integer for a nice y label.
buffer = 3;
uplim = max(uplimrow) + buffer;
lolim = min(lolimrow) - buffer;
uplim = round(uplim);
lolim = round(lolim);

% Create box plot 
% Include significance notches on boxes
% supress outliers
figure(1)
boxplot(S_err_plot,all_labels,'notch','on','symbol','')
ylim([lolim,uplim]);
set(gca,'YTick',lolim:0.5:uplim,'XTickLabelRotation',-30);
ytickformat('percentage');
ylabel('Error');
hold on
yline(0);
hold off





% ================ Area Inertia Rind Only ====================

% Varried ellipse data for RIND ONLY
VE = Ir(:,3:numNEPCs+3,:);

% Get percent error for area moment of inertia for each varried ellipse 
% compared to real shape RIND ONLY
for i = 1:numNEPCs+1
    S_err_rind(:,i,:) =  (VE(:,i,:) - Ir(:,2,:))./Ir(:,2,:)*100;
    S_err_rind_prctile_i(i,:,:) = prctile(S_err_rind(:,i,:),[5,25,50,75,95]);
    for j = 1:5
       S_err_rind_prctile(i,j) = mean(S_err_rind_prctile_i(i,:,j)); 
    end
end



% Create labels according to the number of principal components used in
% the study (cumulative cases followed by remaining individual cases) RIND
% ONLY
all_labels_rind = strings(1,(1+numNEPCs));
all_labels_rind(1,1:2) = ["Ellipse","Ellipse + PC 1"];
for i = 2:numNEPCs
    addlabel_rind = "Ellipse + PCs 1-" + num2str(i);
    all_labels_rind(1,i+1) = addlabel_rind;
end

% Combine all percent error data into one matrix for box plotting RIND ONLY
S_err_rind_plot = S_err_rind(:,:,1);
for i = 1:length(slices)-1
S_err_rind_plot = [S_err_rind_plot;S_err_rind(:,:,i+1)];
end

% Grab upper and lower limits by using 5th and 95th percentile data RIND
% ONLY
uplimrow_rind = S_err_rind_prctile(:,5)';
lolimrow_rind = S_err_rind_prctile(:,1)';

% Add a buffer between the calculated outer reach of the whiskers and the
% edge of the plot. Round to the nearest integer for a nice y label. RIND
% ONLY
buffer_rind = 3;
uplim_rind = max(uplimrow_rind) + buffer_rind;
lolim_rind = min(lolimrow_rind) - buffer_rind;
uplim_rind = round(uplim_rind);
lolim_rind = round(lolim_rind);

% Create box plot for RIND ONLY
% Include significance notches on boxes
% supress outliers
figure(2)
boxplot(S_err_rind_plot,all_labels_rind,'notch','on','symbol','')
ylim([lolim_rind,uplim_rind]);
set(gca,'YTick',lolim_rind:0.5:uplim_rind,'XTickLabelRotation',-30);
ytickformat('percentage');
ylabel('Error');
hold on
yline(0);
hold off








% ============== Max Stress (ms) Percent Error ==============


% Varried ellipse data for MAX STRESS
VE_ms = sigma(:,3:numNEPCs+3,:);

% Get percent error for MAX STRESS for each varried ellipse 
% compared to real shape 
for i = 1:numNEPCs+1
    S_err_ms(:,i,:) =  (VE_ms(:,i,:) - sigma(:,2,:))./sigma(:,2,:)*100;
    S_err_ms_prctile_i(i,:,:) = prctile(S_err_ms(:,i,:),[5,25,50,75,95]);
    for j = 1:5
       S_err_ms_prctile(i,j) = mean(S_err_ms_prctile_i(i,:,j)); 
    end
end



% Create labels according to the number of principal components used in
% the study (cumulative cases followed by remaining individual cases) MAX
% STRESS
all_labels_ms = strings(1,(1+numNEPCs));
all_labels_ms(1,1:2) = ["Ellipse","Ellipse + PC 1"];
for i = 2:numNEPCs
    addlabel_ms = "Ellipse + PCs 1-" + num2str(i);
    all_labels_ms(1,i+1) = addlabel_ms;
end

% Combine all percent error data into one matrix for box plotting MAX
% STRESS
S_err_ms_plot = S_err_ms(:,:,1);
for i = 1:length(slices)-1
S_err_ms_plot = [S_err_ms_plot;S_err_ms(:,:,i+1)];
end

% Grab upper and lower limits by using 5th and 95th percentile data MAX
% STRESS
uplimrow_ms = S_err_ms_prctile(:,5)';
lolimrow_ms = S_err_ms_prctile(:,1)';

% Add a buffer between the calculated outer reach of the whiskers and the
% edge of the plot. Round to the nearest integer for a nice y label. MAX
% STRESS
buffer_ms = 7;
uplim_ms = max(uplimrow_ms) + buffer_ms;
lolim_ms = min(lolimrow_ms) - buffer_ms;
uplim_ms = round(uplim_ms);
lolim_ms = round(lolim_ms);

% Create box plot for MAX STRESS
% Include significance notches on boxes
% supress outliers
figure(3)
boxplot(S_err_ms_plot,all_labels_ms,'notch','on','symbol','')
ylim([lolim_ms,uplim_ms]);
set(gca,'YTick',lolim_ms:0.5:uplim_ms,'XTickLabelRotation',-30);
ytickformat('percentage');
ylabel('Error');
title('Maximum Bending Stress')
hold on
yline(0);
hold off

end
