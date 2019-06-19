function transverse_wrapper(range,slicedist)
% FILENAME: transverse_wrapper.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: Wrap the majority of the data production process into a single
% script
% 
% 
% INPUTS:
%       range - A 1 x 2 vector of integers, indicating the starting and
%       ending stalk numbers that will be used (values must be between 1
%       and 990)
%       slicedist - A number value (can be integer or decimal) that
%       signifies the distance in millimeters from the nearest node on the
%       stalk at the chosen slice.  
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
%       - The closest slice to slicedist for a given stalk will be chosen,
%       so when verifying results make sure to look at the distribution of
%       chosen slice positions when choosing extremely large magnitude
%       numbers.
%       - Don't go above a magnitude of 90 at all, since there are
%       very few stalks that have slices that far from the node. A value of
%       0 for slicedist will give the closest slice to the node.
%       - The PCA functions look at the first five NEPCs, so the range
%       should always cover at least 5 stalks (realistically, the range
%       should cover much more)
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

%% Process
% Load in the 360-point DCR version of Jared's table
load StalksDCR_360pts.mat
hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

r1 = num2str(range(1));
r2 = num2str(range(2));
dist_int = num2str(round(abs(slicedist),0));
deci = abs(slicedist) - round(abs(slicedist),0);
dist_deci = num2str(deci);
dist_deci = erase(dist_deci,'0.');

if slicedist > 0
    if strcmp(dist_deci,'0')
        slicepos = strcat('_Above_',dist_int);
    else
        slicepos = strcat('_Above_',dist_int,'_',dist_deci);
    end
elseif slicedist < 0
    if strcmp(dist_deci,'0')
        slicepos = strcat('_Below_',dist_int);
    else
        slicepos = strcat('_Below_',dist_int,'_',dist_deci);
    end
else
    slicepos = strcat('_At_Node');
end

output_prefix = strcat('Stalks_',r1,'_',r2,slicepos);

% Choose the cross-section samples
ChooseSectionsName = strcat(output_prefix,'_Sampled.mat');
ChooseSections('samedist',range,slicedist,Stalk_TableDCR,npoints,ChooseSectionsName)

% Manually find the cross-sections that need to be flipped 180 degrees
FlipName = strcat(output_prefix,'_flip_sections.mat');
find_flip_notches(ChooseSectionsName,FlipName)

while 1
    fixes_needed = input('Does the flip vector need manual correction? Y/N ','s');
    switch fixes_needed
        case 'Y'
            load(FlipName);
            openvar('flip_sections');
            disp('Giving control to keyboard for manual editing of flip variable.');
            disp('Use dbcont command to exit keyboard editing mode.');
            keyboard;
            break
        case 'N'
            break
            
        otherwise
            disp('Not a recognized answer. Please try again.');
    end
end


% Flip the cross-sections that need to be flipped, according to the vector
% of flip indicators
% load(FlipName);
FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
flip_notches(FlipName,ChooseSectionsName,FlippedOutputName);

% Calculate the ellipse fits
EllipseName = strcat(output_prefix,'_Ellipses.mat');
ellipse_fitting_V2(FlippedOutputName,EllipseName);

% Plot the ellipse fits and see if any of them have problems
load(EllipseName);

problem_indices = [];
for i = 1:(range(2) - range(1) + 1)
    i
    polarplot(ELLIPSE_T(i,:),ELLIPSE_R_ext(i,:));
    hold on
    polarplot(ELLIPSE_T(i,:),ELLIPSE_R_int(i,:));
    hold off
    s = input('Enter 1 if the ellipse fit has a problem: ');
    s
    if s == 1
        problem_indices = [problem_indices, i];
    else
        continue
    end    
end

while 1
    fixes_needed = input('Does the ellipse problems vector need manual correction? Y/N ','s');
    switch fixes_needed
        case 'Y'
            load(FlipName);
            openvar('problem_indices');
            disp('Giving control to keyboard for manual editing of flip variable.');
            disp('Use dbcont command to exit keyboard editing mode.');
            keyboard;
            break
        case 'N'
            break
            
        otherwise
            disp('Not a recognized answer. Please try again.');
    end
end

problem_indices

NEPCName = strcat(output_prefix,'_PCA.mat');

if isempty(problem_indices)
    % Run PCA
    PCA_ellipse_fits(EllipseName,NEPCName);
    
    % Create the Abaqus Python scripts
    create_cases(NEPCName,EllipseName,ChooseSectionsName,problem_indices,5);
else
    % Remove the problem ellipses and then run PCA again
    GoodEllipseFits = strcat(output_prefix,'_GoodEllipses.mat');
    remove_problem_ellipses(EllipseName,problem_indices,GoodEllipseFits);
    
    % Run PCA
    PCA_ellipse_fits(GoodEllipseFits,NEPCName);
    
    % Create the Abaqus Python scripts
    create_cases(NEPCName,GoodEllipseFits,ChooseSectionsName,problem_indices,5);
end

set(0,'DefaultFigureWindowStyle','normal');

end



%% Localizing all functions used
function ChooseSections(method,range,dist,Table,npoints,SaveName)
% ChooseSections.m: Determine the cross-sections to compile, which is
% determined by a method

% range: For Stalk_Table, this must be a row vector of two integer values
% from 1 to 990.

allrows = size(Table,1);

switch method
    % Choose a number of cross-sections that are all at the same distance
    % from the node
    case 'samedist'
        indices = zeros(1,(range(2) - range(1) + 1));
%         dist = input('Choose approximate slice distance to use: ');
        
        % Step through stalks of interest (defined by range values)
        stalk = 1;
        for i = range(1):range(2)
            A = (Table.StkNum == i);
            
            % Get index of first row that is part of the current stalk
            for j = 1:length(A)
                if A(1) == 1
                    idx_first = 1;
                    break
                elseif A(j) == 1 && A(j-1) == 0
                    idx_first = j;
                end
            end
            
            % Get index of last row that is part of the current stalk
            for j = 2:length(A)
                if A(j) == 0 && A(j-1) == 1
                    idx_last = j-1;
                    break
                end
            end
            
            % Iterate through the indices that are part of the current
            % stalk and determine the index that is closest to the desired
            % value
            diffs = NaN(size(Table.StkNum));            
            for j = idx_first:idx_last
                diffs(j) = dist - Table.SlP(j);
            end
            
            % Get index of closest slice
            [~,closestIndex] = min(abs(diffs));
            
            % Add the index to indices output
            indices(stalk) = closestIndex;
            stalk = stalk + 1;
        end
        
        % Create table from indices for later reference
        selectedTable = Table(indices,:);
        
        % Save compiled slices in arrays for downstream use
        ext_X =     makearray(selectedTable,'Ext_X',npoints);
        ext_Y =     makearray(selectedTable,'Ext_Y',npoints);
        int_X =     makearray(selectedTable,'Int_X',npoints);
        int_Y =     makearray(selectedTable,'Int_Y',npoints);
        ext_T =     makearray(selectedTable,'Ext_T',npoints);
        ext_Rho =   makearray(selectedTable,'Ext_Rho',npoints);
        int_T =     makearray(selectedTable,'Int_T',npoints);
        int_Rho =   makearray(selectedTable,'Int_Rho',npoints);
        avg_rind_thick = Table.rind_t(range(1):range(2));
        
        
        % Output all variables into mat file
        FolderName = pwd;
        SaveFile = fullfile(FolderName, SaveName);
        save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_T','ext_Rho',...
            'int_T','int_Rho','avg_rind_thick','indices','selectedTable','npoints');
        
    case 'wholestalk'
        % Choose a range of stalk numbers, and all the slices from each of
        % the chosen stalks will be chosen
        
    otherwise
        disp('Unknown method.');
end


end


function [exported_array] = makearray(Table,Variable,npoints)
    % Make an array of the chosen variable in the format required to work
    % with the downstream process in ellipse_fittingV1.m and
    % PCA_ellipse_fits.m
    num_slices = size(Table,1);
    exported_array = NaN(npoints,num_slices);
    
    if iscell(Table.(Variable))
        for i = 1:num_slices
            exported_array(:,i) = cell2mat(Table.(Variable)(i));
        end
    else
        for i = 1:num_slices
            exported_array(i,:) = Table.(Variable)(i);
        end
    end
        

end


function find_flip_notches(ChooseSectionsOutput,SaveName)
% FILENAME: find_flip_notches.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: Manually identify the cross-sections that need to be flipped 180
% degrees so the notch is on the left side. Cycles through the data
% contained in a .mat file output from ChooseSections.m, and allows the
% user to tag cross-sections that need to be flipped by entering "1" and
% then Enter. All other sections can be skipped through by hitting Enter.
% 
% 
% INPUTS:
%       ChooseSectionsOutput - A .mat file produced by ChooseSections.m
%       
% OUTPUTS:
%       flip_sections - A vector of 1s and empties (or non-1s) that must be
%       fed into flip_notches.m.
%
%
% NOTES: The reason find_flip_notches.m and flip_notches.m are not combined
% into a single function is so the user has the opportunity to correct for
% miskeyed cross-sections. flip_sections can be manually edited before
% being fed into flip_notches in case the user double-typed or had some
% other problem, since find_flip_notches.m only marches forward through the
% cross-sections, with no opportunity to make corrections during the
% process.
% 
% 
% VERSION HISTORY:
% V1 - Made into a function that works with the updated process flow
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load(ChooseSectionsOutput);

R_ext = ext_Rho;
R_int = int_Rho;
T = ext_T(:,1)';

N = size(R_ext,2);
flip_sections = zeros(N,1);

for i = 1:N
    % Plot each cross section to see if it needs to be flipped 180 degrees
    polarplot(T,R_ext(:,i));
    i
    s = input('Enter 1 if cross section needs to flip: ');
    if isempty(s)
        s = 0;
    end
    flip_sections(i) = s;
%     pause(); 
end


% Save data as mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, SaveName);
save(SaveFile,'flip_sections');


end


function flip_notches(Flip_Indices,ChooseSectionsOutput,FlippedOutputName)
% FILENAME: find_flip_notches.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: 
% 
% 
% INPUTS:
%       flip_sections - The vector of 1s and empties that identifies the
%       cross-sections needing to be flipped 180 degrees
%
%       ChooseSectionsOutput - A .mat file produced by ChooseSections.m
%
%       SaveName - The name for the output .mat file. Make sure to end the
%       name with FLIPPED for consistency.
%       
% OUTPUTS:
%       
%
% NOTES:
%       Make sure to name the output .mat file with a FLIPPED extension
%       for consistency
% 
% 
% VERSION HISTORY:
% V1 - Made into a function that works with the updated process flow
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load(ChooseSectionsOutput);
load(Flip_Indices);

ext_rho = [];
int_rho = [];
% ext_X = ext_X;
% ext_Y = ext_Y;
% int_X = int_X;
% int_Y = int_Y;
T = ext_T(:,1)';
N = size(ext_X,2);

flipped = []; % Hold the indices of the flipped sections
for i = 1:N
    if flip_sections(i) == 1 
        flipped = [flipped; i];
        
        % Rotate and reorder external and internal points
        [~, ~, ~, ~, ~, ~, xp_ext, yp_ext, ~, ~] = reorder_V2(ext_X(:,i), ext_Y(:,i), pi);
        [~, ~, ~, ~, ~, ~, xp_int, yp_int, ~, ~] = reorder_V2(int_X(:,i), int_Y(:,i), pi);
        
        [~, ~, x_ext, y_ext, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_ext, yp_ext, 0);
        [~, ~, x_int, y_int, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_int, yp_int, 0);
        
        % Redefine the appropriate row in the main XY data
        ext_X(:,i) = x_ext;
        ext_Y(:,i) = y_ext;
        int_X(:,i) = x_int;
        int_Y(:,i) = y_int;
    end
end

% % Plot the cross sections that should be flipped to verify that the notches
% % are all on the left side
% for i = 1:length(flipped)
%     for j = 1:length(T)
%         plot(ext_X(1:j,flipped(i)),ext_Y(1:j,flipped(i)));
%         hold on
%         plot(int_X(1:j,flipped(i)),int_Y(1:j,flipped(i)));
%         axis equal
%         hold off
% %         pause(0.01);
%     end
%     pause();
% end


% Convert data from Cartesian to polar
for i = 1:N
    for j = 1:length(T)
        ext_Rho(i,j) = sqrt(ext_X(j,i)^2 + ext_Y(j,i)^2);
        int_Rho(i,j) = sqrt(int_X(j,i)^2 + int_Y(j,i)^2);
    end
end

% Transpose rho arrays so they are the same orientation as the other
% variables
ext_rho = ext_rho';
int_rho = int_rho';

% Save data as mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, FlippedOutputName);
save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_Rho','int_Rho','ext_T','int_T','avg_rind_thick','flip_sections','npoints');

end


function ellipse_fitting_V2(FileName,SaveName)
% Load a mat file that has exterior XY data and avgrindthickness data.
% Cycle through cross sections and select the angular range that
% contains the notch so it's ignored during ellipse fitting. Then save the
% ellipses in a mat file. Take the difference between the interior and
% exterior points and their ellipse approximations and save those as well
% for later PCA.

close all

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'avg_rind_thick','ext_Rho','ext_T','ext_X','ext_Y',...
    'int_Rho','int_T','int_X','int_Y','npoints');

% Make copies of original data to work with
R_ext = ext_Rho;
R_int = int_Rho;
T = (ext_T(:,1));
X_ext = ext_X;
Y_ext = ext_Y;
X_int = int_X;
Y_int = int_Y;

N = length(avg_rind_thick);

% Instantiate output variables
A = zeros(N,1);
B = zeros(N,1);
ALPHA = zeros(N,1);

ELLIPSE_XY = zeros(N,npoints,2);
ELLIPSE_CENTERS = zeros(N,2);
ELLIPSE_T = zeros(N,npoints);
ELLIPSE_R_ext = zeros(N,npoints);
ELLIPSE_R_int = zeros(N,npoints);

DIFF_R_ext = zeros(N,npoints);
DIFF_R_int = zeros(N,npoints);
R_ext = zeros(N,npoints);
R_int = zeros(N,npoints);

AVG_RIND_T = zeros(N,1);

prev_alpha = 0;

min_angle = 135;
min_angle = min_angle*(pi/180);     % Convert angle to radians
max_angle = 225;
max_angle = max_angle*(pi/180);

for i = 1:N
    % Define the notch range
    for j = 1:npoints
        if T(j) > min_angle
            min_index = j-1;
            break
        end
    end
    
    for j = 1:npoints
        if T(j) > max_angle
            max_index = j;
            break
        end
    end
    
    % Cut out the notch from the XY data
    window = [linspace(1,min_index,min_index),linspace(max_index,npoints,(npoints-max_index + 1))];
    xcut = X_ext(window,i);
    ycut = Y_ext(window,i);

    % Fit an ellipse to the data with the notch removed
    [alpha, major, minor, xbar_e, ybar_e, X_ellipse, Y_ellipse] = fit_ellipse_R4( xcut, ycut, npoints, prev_alpha, gca );
    
    % Save ellipse center shift
    ELLIPSE_CENTERS(i,:) = [mean(X_ellipse), mean(Y_ellipse)];
    
    % Shift XY data to be centered at the origin before converting to polar
    X_ellipse_shift = X_ellipse - mean(X_ellipse);
    Y_ellipse_shift = Y_ellipse - mean(Y_ellipse);
    X_ext_shift = X_ext(:,i) - mean(X_ellipse);
    Y_ext_shift = Y_ext(:,i) - mean(Y_ellipse);
    X_int_shift = X_int(:,i) - mean(X_ellipse);
    Y_int_shift = Y_int(:,i) - mean(Y_ellipse);
    
    % Reorder indices to start at 0 degrees
    [X_ellipse_shift, Y_ellipse_shift, ~, ~, ~, ~, ~, ~, ~, ~] = reorder_V2(X_ellipse_shift, Y_ellipse_shift, 0);
    [X_ext_shift, Y_ext_shift, ~, ~, ~, ~, ~, ~, ~, ~] = reorder_V2(X_ext_shift, Y_ext_shift, 0);
    [X_int_shift, Y_int_shift, ~, ~, ~, ~, ~, ~, ~, ~] = reorder_V2(X_int_shift, Y_int_shift, 0);    
    
    % Convert X_ellipse and Y_ellipse to polar coordinates
    theta = 0:2*pi/npoints:2*pi;
    theta = theta(1:end-1);
    [thetatemp_ellipse, ext_rho_ellipse] = cart2pol(X_ellipse_shift,Y_ellipse_shift);
    thetatemp_ellipse = wrapTo2Pi(thetatemp_ellipse);   % Make all the negative pi values sit on a 0-2*pi system
    
    % Remove duplicate values for interpolation
    [C,ia,~] = unique(thetatemp_ellipse,'stable');
    thetatemp_ellipse = C;
    ext_rho_ellipse = ext_rho_ellipse(ia);
    
    % Convert X_ext_shift and Y_ext_shift to polar for resampling
    [thetatemp_ext, ext_rho] = cart2pol(X_ext_shift,Y_ext_shift);
    thetatemp_ext = wrapTo2Pi(thetatemp_ext);   % Make all the negative pi values sit on a 0-2*pi system
    
    % Remove duplicate values for interpolation
    [C,ia,~] = unique(thetatemp_ext,'stable');
    thetatemp_ext = C;
    ext_rho = ext_rho(ia);
    
    % Convert X_int_shift and Y_int_shift to polar for resampling
    [thetatemp_int, int_rho] = cart2pol(X_int_shift,Y_int_shift);
    thetatemp_int = wrapTo2Pi(thetatemp_int);   % Make all the negative pi values sit on a 0-2*pi system
    
    % Remove duplicate values for interpolation
    [C,ia,~] = unique(thetatemp_int,'stable');
    thetatemp_int = C;
    int_rho = int_rho(ia);
    
    
    
    
    % Interpolate to get new rho and theta points that are regularly spaced
    ext_rho_ellipse_interp = interp1(thetatemp_ellipse,ext_rho_ellipse,theta,'pchip','extrap'); 
    ext_rho_ellipse = ext_rho_ellipse_interp;
    ext_rho_interp = interp1(thetatemp_ext,ext_rho,theta,'pchip','extrap'); 
    ext_rho = ext_rho_interp;
    int_rho_interp = interp1(thetatemp_int,int_rho,theta,'pchip','extrap'); 
    int_rho = int_rho_interp;
    
    % Get interior ellipse fit points, based on constant rind thickness
    % assumption
    int_rho_ellipse = ext_rho_ellipse - avg_rind_thick(i);
    
%     % Plot in polar coordinates to check results
%     polarplot(theta,ext_rho,'.','LineWidth',2);
%     hold on
%     polarplot(theta,int_rho,'.','LineWidth',2);
% %     polarplot(theta,ext_rho_ellipse,'.','LineWidth',2)
%     pause();
%     close;
    
    
    
    A(i) = major;
    B(i) = minor;
    ALPHA(i) = alpha;

    ELLIPSE_XY(i,:,1) = X_ellipse;
    ELLIPSE_XY(i,:,2) = Y_ellipse;
    ELLIPSE_CENTERS(i,1) = mean(X_ellipse);
    ELLIPSE_CENTERS(i,2) = mean(Y_ellipse);
    ELLIPSE_T(i,:) = theta;
    ELLIPSE_R_ext(i,:) = ext_rho_ellipse;
%     ELLIPSE_R_int(i,:) = ext_rho_ellipse - avg_rind_thick(i);
    ELLIPSE_R_int(i,:) = rpts(npoints,ELLIPSE_T(i,:),(A(i) - 2*avg_rind_thick(i)),(B(i) - 2*avg_rind_thick(i)));
    R_ext(i,:) = ext_rho;
    R_int(i,:) = int_rho;
    
    % Get difference between the ellipse and the real data (if the ellipse
    % overestimates, then the value of DIFF will be positive)
    DIFF_R_ext(i,:) = ext_rho_ellipse - ext_rho;
    DIFF_R_int(i,:) = int_rho_ellipse - int_rho;
    
    AVG_RIND_T(i) = avg_rind_thick(i);
    
    RIND_ELLIPSE_DIFF(i,:) = ELLIPSE_R_ext(i,:) - R_int(i,:);
    
end

% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int',...
    'ELLIPSE_CENTERS','DIFF_R_ext','DIFF_R_int','R_ext','R_int','AVG_RIND_T');

end


function [r] = rpts(N,theta,dmaj,dmin)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end


function PCA_ellipse_fits(FileName,SaveName)
% USE THIS FUNCTION ON Ellipse_fits_bottom1.mat or Ellipse_fits_top1.mat
% (uses the difference between the ellipse and the real data)

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'DIFF_R_ext','DIFF_R_int','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int');

% Perform PCA. 'Centered' option must be false to allow for reverse
% engineering of the original data
[ext_rhoPCAs, ext_rhocoeffs, ext_rhoPCA_variances, ext_rhotstat, ext_rhoexplained, ext_rhovarMeans] = pca(DIFF_R_ext,'Centered',false);
% [int_rhoPCAs, int_rhocoeffs, int_rhoPCA_variances, int_rhotstat, int_rhoexplained, int_rhovarMeans] = pca(DIFF_R_int,'Centered',false);

ext_rhoexplained_tot = zeros(size(ext_rhoexplained));
% int_rhoexplained_tot = zeros(size(ext_rhoexplained));
for i = 1:length(ext_rhoexplained_tot)
    ext_rhoexplained_tot(i) = sum(ext_rhoexplained(1:i));
%     int_rhoexplained_tot(i) = sum(int_rhoexplained(1:i));
end

figure(1);
plot(ext_rhoexplained_tot(:,1),'-*');
title('Exterior Non-Ellipse PCs');
xlabel('# of PCs');
ylabel('% Variance Explained');

% figure(2);
% plot(int_rhoexplained_tot(:,1),'-*');
% title('Interior Non-Ellipse PCs');
% xlabel('# of PCs');
% ylabel('% Variance Explained');

ELLIPSE_T = ELLIPSE_T';
theta = ELLIPSE_T(:,1);
ELLIPSE_T = ELLIPSE_T';

figure(3);
polarplot(theta,ext_rhoPCAs(:,1));
hold on
polarplot(theta,ext_rhoPCAs(:,2));
polarplot(theta,ext_rhoPCAs(:,3));
polarplot(theta,ext_rhoPCAs(:,4));
polarplot(theta,ext_rhoPCAs(:,5));
title('Exterior Rho Principal Components');
legend('PC1','PC2','PC3','PC4','PC5');

% figure(4);
% polarplot(theta,int_rhoPCAs(:,1));
% hold on
% polarplot(theta,int_rhoPCAs(:,2));
% polarplot(theta,int_rhoPCAs(:,3));
% polarplot(theta,int_rhoPCAs(:,4));
% polarplot(theta,int_rhoPCAs(:,5));
% title('Interior Rho Principal Components');
% legend('PC1','PC2','PC3','PC4','PC5');




% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'ELLIPSE_T','ELLIPSE_R_ext','ext_rhocoeffs',...
    'ext_rhoPCAs','ext_rhoexplained','ext_rhovarMeans');


end


function remove_problem_ellipses(OriginalEllipseFits,problem_indices,GoodEllipseFits)
load(OriginalEllipseFits);

N = size(A,1);

% Remove the problem ellipses from the ellipse fit data
A(problem_indices) = [];
AVG_RIND_T(problem_indices) = [];
B(problem_indices) = [];
DIFF_R_ext(problem_indices,:) = [];
DIFF_R_int(problem_indices,:) = [];
ELLIPSE_CENTERS(problem_indices,:) = [];
ELLIPSE_R_ext(problem_indices,:) = [];
ELLIPSE_R_int(problem_indices,:) = [];
ELLIPSE_T(problem_indices,:) = [];
ELLIPSE_XY(problem_indices,:,:) = [];
R_ext(problem_indices,:) = [];
R_int(problem_indices,:) = [];

% Save the final data in a new mat file
FolderName = pwd;
SaveFile       = fullfile(FolderName, GoodEllipseFits);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int',...
    'ELLIPSE_CENTERS','DIFF_R_ext','DIFF_R_int','R_ext','R_int','AVG_RIND_T','problem_indices');

end



