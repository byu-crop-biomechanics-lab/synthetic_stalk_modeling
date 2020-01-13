function PrepSections3D(stalknums, npoints, DataTable, SaveName)
% FILENAME: PrepSections3D.m
% AUTHOR: Ryan Larson
% DATE: 1/9/2020
%
% PURPOSE: 
%   This is the parallel to PrepSections_V2.m for the transverse model
%   process. This function takes in the table prepared by Jared that has
%   exterior and interior stalk cross-section boundaries and shifts each
%   stalk so the node cross-section is centered and rotated. Cross-sections
%   below the node will have the notch on the left, and sections above the
%   node will have the notch on the right. The data is downsampled to
%   include 360 points.
%   THE TABLE WILL ALSO REMOVE ANY DATA POINTS THAT ARE FOUND TO CAUSE A
%   PROCESSING ERROR, AS CAUGHT BY THE TRY CATCH BLOCK. THIS IS TO PREVENT
%   PROBLEMS DOWN THE LINE WITH PROCESSING.
% 
% INPUTS: 
%       stalknums - A 2x1 vector defining the starting and ending stalk
%       numbers to be included in the output data table
%
%       npoints - An integer value for the number of evenly-spaced points
%       in polar coordinates to use. This does not include a repeat of the
%       first point to close the shape, so if you want the points to be
%       exactly on 1-degree intervals, set npoints = 359.
%
%       DataTable - This should be the data table Stalk_Table from
%       SMALL_CURVES_V2_3_1500.mat
%
%       SaveName - A string with an output file name, with .mat file
%       extension.
%
%
% OUTPUTS:
%   This function outputs a .mat file containing the following variables:
%       Stalk_TableDCR - A table that is the same as Stalk_Table from the
%       input, except that Ext_X, Ext_Y, Int_X, and Int_Y are replaced with
%       their downsampled, centered, and rotated versions. Added are Ext_T,
%       Ext_Rho, Int_T, and Int_Rho, the polar coordinate versions of the 
%       downsampled data.
% 
%       error_indices - A list of integers corresponding to rows in
%       Stalk_Table that experienced errors during conversion, and should
%       not be used for further analysis until the problems are fixed.
%       ACTUALLY REMOVE THESE INDICES BEFORE SAVING THE TABLE BECAUSE
%       THEY'RE A
%       HEADACHE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%       npoints - The same integer value from the input to this function,
%       saved to ensure that all downstream functions assume the same
%       number of downsampled cross-section points.
%
%
% NOTES: - Originally adapted from PrepSections_V2.m
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

% Initialize variables
plotting = 0;       % a following function has a built-in plotting option, which we turn off

%%% Variable Initializations
alpha = 0;
prev_alpha = 0;

nslices = 0; % Initialize the total number of slices in the chosen stalks
startslice = 1;
endslice = 1;

if stalknums(1) ~= 1
    for i = 1:size(DataTable,1)-1
        if DataTable.StkNum(i) == stalknums(1) && DataTable.StkNum(i-1) ~= DataTable.StkNum(i)
            startslice = i;
        end
    end
else
    startslice = 1;
end

for i = 1:size(DataTable,1)-1
    if DataTable.StkNum(i) == stalknums(2) && DataTable.StkNum(i+1) ~= DataTable.StkNum(i)
        endslice = i;
    end
end

nslices = endslice - startslice + 1;


error_indices = []; % Global indices where errors occur
theta_rot = zeros(nslices,1);
A = zeros(nslices,1);
B = zeros(nslices,1);

nstalks = stalknums(2) - stalknums(1) + 1;
nodeindices = zeros(nstalks,1); % Holding vector for global indices of nodes

Stalk_TablePCA = DataTable(stalknums(1):stalknums(2),:);

%% Start loop
% Loop through stalk numbers
for n = stalknums(1):stalknums(2)
    n
    % Get the starting and ending indices of the table for the current stalk
    % number
    indices = (DataTable.StkNum == n); % This outputs a logical where the condition is true, not a table
            
    % Get index of first row that is part of the current stalk
    for i = 1:length(indices)
        if indices(1) == 1
            idx_first = 1;
            break
        elseif indices(i) == 1 && indices(i-1) == 0
            idx_first = i;
        end
    end

    % Get index of last row that is part of the current stalk
    for i = 2:length(indices)
        if indices(i) == 0 && indices(i-1) == 1
            idx_last = i-1;
            break
        end
    end
    
    tempTable = DataTable((idx_first:idx_last),:);

    % Locate node cross-section in the complete original DataTable, not the
    % indices selection
    diffs = NaN(size(DataTable.StkNum));            
    for i = idx_first:idx_last
        diffs(i) = 0 - DataTable.SlP(i);
    end
    
    % Get index of node on full table level
    [~,globalNodeIndex] = min(abs(diffs));
    nodeindices(n) = globalNodeIndex; % Save the global node indices
    
    % Locate node cross-section in the sub-table for just this stalk
    diffs = NaN(size(tempTable.StkNum));            
    for i = 1:size(tempTable,1)
        diffs(i) = 0 - tempTable.SlP(i);
    end
    
    % Get index of node on current stalk table level
    [~,localNodeIndex] = min(abs(diffs));
    nodeindices(n) = localNodeIndex; % Save the global node indices
    
    
    % Define geometric center of the node cross-section
    xcnode = DataTable.xbar(globalNodeIndex);
    ycnode = DataTable.ybar(globalNodeIndex);
    
    % Shift all cross-sections in the current stalk by the node
    % cross-section shift
    xshift = xcnode;
    yshift = ycnode;
    for i = 1:size(tempTable,1)
        % Convert cells to arrays for easier access
        ext_X = cell2mat(tempTable.Ext_X(i));
        ext_Y = cell2mat(tempTable.Ext_Y(i));
        int_X = cell2mat(tempTable.Int_X(i));
        int_Y = cell2mat(tempTable.Int_Y(i));
        
        % Shift current cross-section by xshift and yshift
        for j = 1:length(ext_X)
            ext_X(j) = ext_X(j) - xshift;
        end
        for j = 1:length(ext_Y)
            ext_Y(j) = ext_Y(j) - yshift;
        end
        for j = 1:length(int_X)
            int_X(j) = int_X(j) - xshift;
        end
        for j = 1:length(int_Y)
            int_Y(j) = int_Y(j) - yshift;
        end
        
        tempTable.Ext_X(i) = {ext_X};
        tempTable.Ext_Y(i) = {ext_Y};
        tempTable.Int_X(i) = {int_X};
        tempTable.Int_Y(i) = {int_Y};
        
        tempTable.xbar(i) = mean(ext_X);
        tempTable.ybar(i) = mean(ext_Y);
        
    end
    
    %% At this point, a smaller table exists for the current stalk only, with all cross-sections shifted together so that the node cross-section is centered at the origin
    % All operations should be done on this, then inserted into the
    % original table at the appropriate rows.
    
    %% Rotate all cross-sections about the center of the node cross-section  (which should be at 0,0)
    
    % Determine the angle of the node cross-section
    prev_alpha = 0;
    ext_X = cell2mat(tempTable.Ext_X(localNodeIndex));
    ext_Y = cell2mat(tempTable.Ext_Y(localNodeIndex));
    
    [node_alpha,~,~] = getrotation(ext_X, ext_Y, prev_alpha);
    
    % Check the notch angle of the first 10 cross-sections. Determine the
    % correct hemisphere to rotate the node into by taking the rough angle
    % of the majority of these and determining what quadrant will put the
    % notch at the left for these cross-sections. Then use this as a
    % correction if the notch is difficult to locate for the node. This
    % should avoid whole stalks being turned the wrong way.
    nbottom = 4;
    bottom_angles = zeros(nbottom,1);
    notch_indicator = zeros(nbottom,1);
    for i = 1:nbottom
        prev_alpha = 0;
        ext_X = cell2mat(tempTable.Ext_X(i));
        ext_Y = cell2mat(tempTable.Ext_Y(i));
        [bottom_angles(i),~,~] = getrotation(ext_X, ext_Y, prev_alpha);
        [~, ~, ~, ~, ~, ~, ext_xi, ext_yi, ~, ~] = reorder_V2(ext_X, ext_Y, bottom_angles(i));

        plot(ext_xi,ext_yi);
        axis equal
        
        s = input('Enter 1 if notch is NOT on the left: ');
        if isempty(s) || s ~= 1
            s = 0;
        end
        notch_indicator(i) = s;
    end
    
    
    % Add rotation to node_alpha if the notch at the bottom of the stalk is
    % not on the left
    
    count_flip = sum(notch_indicator(:) == 1);
    
    if count_flip > 2
        if node_alpha > 0
            node_alpha = node_alpha + pi;
        else
            node_alpha = node_alpha - pi;
        end
    end
    
    % Rotate all exteriors and interiors by node_alpha, about 0,0
    for i = 1:size(tempTable,1)
        try
        
            ext_X = cell2mat(tempTable.Ext_X(i));
            ext_Y = cell2mat(tempTable.Ext_Y(i));
            int_X = cell2mat(tempTable.Int_X(i));
            int_Y = cell2mat(tempTable.Int_Y(i));
            [~, ~, ~, ~, ~, ~, ext_xi, ext_yi, ~, ~] = reorder_V2_interior(ext_X, ext_Y, node_alpha, 0, 0);
            [~, ~, ~, ~, ~, ~, int_xi, int_yi, ~, ~] = reorder_V2_interior(int_X, int_Y, node_alpha, 0, 0);

            % Save rotation angles of each cross-section relative to the x-axis.
            prev_alpha = 0;
    %         ext_X = cell2mat(tempTable.Ext_X(i));
    %         ext_Y = cell2mat(tempTable.Ext_Y(i));        
            [theta_rot(i),A(i),B(i)] = getrotation(ext_xi, ext_yi, prev_alpha);

            % Downsample cross-sections
            idx =  1:length(ext_xi);                            % Index
            idxq = linspace(min(idx), max(idx), npoints);       % Interpolation Vector
            ext_xi = interp1(idx, ext_xi, idxq, 'pchip');       % Downsampled Vector

            idy = 1:length(ext_yi);                             % Index
            idyq = linspace(min(idy), max(idy), npoints);       % Interpolation Vector
            ext_yi = interp1(idy, ext_yi, idyq, 'pchip');       % Downsampled Vector

            idx =  1:length(int_xi);                            % Index
            idxq = linspace(min(idx), max(idx), npoints);       % Interpolation Vector
            int_xi = interp1(idx, int_xi, idxq, 'pchip');       % Downsampled Vector

            idy = 1:length(int_yi);                             % Index
            idyq = linspace(min(idy), max(idy), npoints);       % Interpolation Vector
            int_yi = interp1(idy, int_yi, idyq, 'pchip');       % Downsampled Vector

            % Temporarily shift boundary data to be centered at the
            % centroid of the exterior data so conversion to polar and
            % interpolation goes better.
            xshift = tempTable.xbar(i);
            yshift = tempTable.ybar(i);
            ext_xi = ext_xi - xshift;
            ext_yi = ext_yi - yshift;
            int_xi = int_xi - xshift;
            int_yi = int_yi - yshift;

            % Get interior and exterior data in polar coordinates
            [~, ~, ~, ~, ~, ~, ~, ~, ext_rho, ext_t] = reorder_V2_interior(ext_xi, ext_yi, 0, 0, 0);
            [~, ~, ~, ~, ~, ~, ~, ~, int_rho, int_t] = reorder_V2_interior(int_xi, int_yi, 0, 0, 0);

            % Interpolate in polar to get first point exactly on the x-axis
            % when converted back to Cartesian
            theta = linspace(0,2*pi,npoints+1); % npoints points from 0 to 2*pi inclusive (puts the theta values right on degrees if npoints = 359)
            theta = transpose(theta(1:end-1)); % Remove the last point so there are npoints points in the end

            [ext_t, index] = unique(ext_t); 
            ext_rho_interp = interp1(ext_t,ext_rho(index),theta,'pchip','extrap');
            [int_t, index] = unique(int_t);
            int_rho_interp = interp1(int_t,int_rho(index),theta,'pchip','extrap');  
            ext_t = theta;
            int_t = theta;
            ext_rho = ext_rho_interp;
            int_rho = int_rho_interp;

            % Convert the resampled polar points back to Cartesian before
            % output
            ext_xi = ext_rho.*cos(ext_t);
            ext_yi = ext_rho.*sin(ext_t);
            int_xi = int_rho.*cos(int_t);
            int_yi = int_rho.*sin(int_t);

            % Shift the boundaries back by their shift values
            ext_xi = ext_xi + xshift;
            ext_yi = ext_yi + yshift;
            int_xi = int_xi + xshift;
            int_yi = int_yi + yshift;

            % Insert rotated and downsampled data back into tempTable
            tempTable.Ext_X(i) = {ext_xi};
            tempTable.Ext_Y(i) = {ext_yi};
            tempTable.Int_X(i) = {int_xi};
            tempTable.Int_Y(i) = {int_yi};
        
        % Catch error cases from try block
        catch
            global_index = idx_first + i - 1;
            error_indices = [error_indices, global_index];
        end
        
    end
    
    % Insert tempTable into copy of DataTable
    Stalk_TablePCA(idx_first:idx_last,:) = tempTable;

end

assignin('base','nodeindices',nodeindices);
assignin('base','A',A);
assignin('base','B',B);
assignin('base','theta_rot',theta_rot);
Stalk_TablePCA = addvars(Stalk_TablePCA,A);
Stalk_TablePCA = addvars(Stalk_TablePCA,B);
Stalk_TablePCA = addvars(Stalk_TablePCA,theta_rot);



%% Remove error cases from the table
% Sort error_indices in descending order
error_indices = sort(error_indices,'descend');
for i = 1:length(error_indices)
    Stalk_TablePCA(error_indices(i),:) = [];
end



% Output all variables into mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, SaveName);
save(SaveFile,'Stalk_TablePCA','error_indices','npoints');


end




function [alpha, major, minor, xbar_e, ybar_e, X_ellipse, Y_ellipse] = fit_ellipse_R2( x, y, prev_alpha, axis_handle )
%
% fit_ellipse - finds the best fit to an ellipse for the given set of points.
%
% Format:   ellipse_t = fit_ellipse( x,y,axis_handle )
%
% Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
%           axis_handle - optional. a handle to an axis, at which the estimated ellipse 
%                         will be drawn along with it's axes
%
% OUTPUT:   alpha: Rotation of the ellipse relative to the Cartesian
%                  coordinate system
%           major: Major diameter
%           minor: Minor diameter
%           xbar_e: Center of the ellipse
%           ybar_e: Center of the ellipse
%           X_ellipse: X data (relative to original coordinates?)
%           Y_ellipse: Y data
%
% Note:     if an ellipse was not detected (but a parabola or hyperbola), then
%           an empty structure is returned
%  
% IMPORTANT NOTE: alpha values are based on IMAGE coordinates, in which x is
% horizontal and the y axis points DOWN! This means that a positive alpha
% is a CW angle from the horizontal to the major axis!!!
%
% =====================================================================================
%                  Ellipse Fit using Least Squares criterion
% =====================================================================================
% We will try to fit the best ellipse to the given measurements. the mathematical
% representation of use will be the CONIC Equation of the Ellipse which is:
% 
%    Ellipse = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
%   
% The fit-estimation method of use is the Least Squares method (without any weights)
% The estimator is extracted from the following equations:
%
%    g(x,y;A) := a*x^2 + b*x*y + c*y^2 + d*x + e*y = f
%
%    where:
%       A   - is the vector of parameters to be estimated (a,b,c,d,e)
%       x,y - is a single measurement
%
% We will define the cost function to be:
%
%   Cost(A) := (g_c(x_c,y_c;A)-f_c)'*(g_c(x_c,y_c;A)-f_c)
%            = (X*A+f_c)'*(X*A+f_c) 
%            = A'*X'*X*A + 2*f_c'*X*A + N*f^2
%
%   where:
%       g_c(x_c,y_c;A) - vector function of ALL the measurements
%                        each element of g_c() is g(x,y;A)
%       X              - a matrix of the form: [x_c.^2, x_c.*y_c, y_c.^2, x_c, y_c ]
%       f_c            - is actually defined as ones(length(f),1)*f
%
% Derivation of the Cost function with respect to the vector of parameters "A" yields:
%
%   A'*X'*X = -f_c'*X = -f*ones(1,length(f_c))*X = -f*sum(X)
%
% Which yields the estimator:
%
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       |  A_least_squares = -f*sum(X)/(X'*X) ->(normalize by -f) = sum(X)/(X'*X)  |
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% (We will normalize the variables by (-f) since "f" is unknown and can be accounted for later on)
%  
% NOW, all that is left to do is to extract the parameters from the Conic Equation.
% We will deal the vector A into the variables: (A,B,C,D,E) and assume F = -1;
%
%    Recall the conic representation of an ellipse:
% 
%       A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
% 
% We will check if the ellipse has a tilt (=orientation). The orientation is present
% if the coefficient of the term "x*y" is not zero. If so, we first need to remove the
% tilt of the ellipse.
%
% If the parameter "B" is not equal to zero, then we have an orientation (tilt) to the ellipse.
% we will remove the tilt of the ellipse so as to remain with a conic representation of an 
% ellipse without a tilt, for which the math is more simple:
%
% Non tilt conic rep.:  A`*x^2 + C`*y^2 + D`*x + E`*y + F` = 0
%
% We will remove the orientation using the following substitution:
%   
%   Replace x with cx+sy and y with -sx+cy such that the conic representation is:
%   
%   A(cx+sy)^2 + B(cx+sy)(-sx+cy) + C(-sx+cy)^2 + D(cx+sy) + E(-sx+cy) + F = 0
%
%   where:      c = cos(phi)    ,   s = sin(phi)
%
%   and simplify...
%
%       x^2(A*c^2 - Bcs + Cs^2) + xy(2A*cs +(c^2-s^2)B -2Ccs) + ...
%           y^2(As^2 + Bcs + Cc^2) + x(Dc-Es) + y(Ds+Ec) + F = 0
%
%   The orientation is easily found by the condition of (B_new=0) which results in:
% 
%   2A*cs +(c^2-s^2)B -2Ccs = 0  ==> phi = 1/2 * atan( b/(c-a) )
%   
%   Now the constants   c=cos(phi)  and  s=sin(phi)  can be found, and from them
%   all the other constants A`,C`,D`,E` can be found.
%
%   A` = A*c^2 - B*c*s + C*s^2                  D` = D*c-E*s
%   B` = 2*A*c*s +(c^2-s^2)*B -2*C*c*s = 0      E` = D*s+E*c 
%   C` = A*s^2 + B*c*s + C*c^2
%
% Next, we want the representation of the non-tilted ellipse to be as:
%
%       Ellipse = ( (X-X0)/a )^2 + ( (Y-Y0)/b )^2 = 1
%
%       where:  (X0,Y0) is the center of the ellipse
%               a,b     are the ellipse "radiuses" (or sub-axis)
%
% Using a square completion method we will define:
%       
%       F`` = -F` + (D`^2)/(4*A`) + (E`^2)/(4*C`)
%
%       Such that:    a`*(X-X0)^2 = A`(X^2 + X*D`/A` + (D`/(2*A`))^2 )
%                     c`*(Y-Y0)^2 = C`(Y^2 + Y*E`/C` + (E`/(2*C`))^2 )
%
%       which yields the transformations:
%       
%           X0  =   -D`/(2*A`)
%           Y0  =   -E`/(2*C`)
%           a   =   sqrt( abs( F``/A` ) )
%           b   =   sqrt( abs( F``/C` ) )
%
% And finally we can define the remaining parameters:
%
%   long_axis   = 2 * max( a,b )
%   short_axis  = 2 * min( a,b )
%   Orientation = phi
%
%

% initialize
orientation_tolerance = 1e-3;

% empty warning stack
warning( '' );

% prepare vectors, must be column vectors
x = x(:);
y = y(:);

% remove bias of the ellipse - to make matrix inversion more accurate. (will be added later on).
mean_x = mean(x);
mean_y = mean(y);
x = x-mean_x;
y = y-mean_y;

% the estimation for the conic equation of the ellipse
X = [x.^2, x.*y, y.^2, x, y ];
a = sum(X)/(X'*X);

% check for warnings
if ~isempty( lastwarn )
    disp( 'stopped because of a warning regarding matrix inversion' );
    ellipse_t = [];
    return
end

% extract parameters from the conic equation
[a,b,c,d,e] = deal( a(1),a(2),a(3),a(4),a(5) );

% remove the orientation from the ellipse
if ( min(abs(b/a),abs(b/c)) > orientation_tolerance )
    
    orientation_rad = 1/2 * atan2( b,(c-a) );
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    [a,b,c,d,e] = deal(...
        a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
        0,...
        a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
        d*cos_phi - e*sin_phi,...
        d*sin_phi + e*cos_phi );
    [mean_x,mean_y] = deal( ...
        cos_phi*mean_x - sin_phi*mean_y,...
        sin_phi*mean_x + cos_phi*mean_y );
else
    orientation_rad = 0;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
end

% check if conic equation represents an ellipse
test = a*c;
switch (1)
case (test>0),  status = '';
case (test==0), status = 'Parabola found';  warning( 'fit_ellipse: Did not locate an ellipse' );
case (test<0),  status = 'Hyperbola found'; warning( 'fit_ellipse: Did not locate an ellipse' );
end

% if we found an ellipse return it's data
if (test>0)
    
    % make sure coefficients are positive as required
    if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end
    
    % final ellipse parameters
    X0          = mean_x - d/2/a;
    Y0          = mean_y - e/2/c;
    F           = 1 + (d^2)/(4*a) + (e^2)/(4*c);
    [a,b]       = deal( sqrt( F/a ),sqrt( F/c ) );    
    long_axis   = 2*max(a,b);
    short_axis  = 2*min(a,b);

    % rotate the axes backwards to find the center point of the original TILTED ellipse
    R           = [ cos_phi sin_phi; -sin_phi cos_phi ];
    P_in        = R * [X0;Y0];
    X0_in       = P_in(1);
    Y0_in       = P_in(2);
    
    % pack ellipse into a structure
    ellipse_t = struct( ...
        'a',a,...
        'b',b,...
        'phi',orientation_rad,...
        'X0',X0,...
        'Y0',Y0,...
        'X0_in',X0_in,...
        'Y0_in',Y0_in,...
        'long_axis',long_axis,...
        'short_axis',short_axis,...
        'status','' );
else
    % report an empty structure
    ellipse_t = struct( ...
        'a',[],...
        'b',[],...
        'phi',[],...
        'X0',[],...
        'Y0',[],...
        'X0_in',[],...
        'Y0_in',[],...
        'long_axis',[],...
        'short_axis',[],...
        'status',status );
end

% check if we need to plot an ellipse with its axes.
if (nargin>2) & ~isempty( axis_handle ) & (test>0) & axis_handle~=0
    
    % rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    
    % the axes
    ver_line        = [ [X0 X0]; Y0+b*[-0.75 0.75] ];
    horz_line       = [ X0+a*[-0.75 0.75]; [Y0 Y0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    
    % the ellipse
    theta_r         = linspace(0,2*pi,360);
    ellipse_x_r     = X0 + a*cos( theta_r );
    ellipse_y_r     = Y0 + b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    % draw
    %hold_state = get( axis_handle,'NextPlot' );
    %set( axis_handle,'NextPlot','add' );
%     plot( new_ver_line(1,:),new_ver_line(2,:),'r' ,'LineWidth', 1);
%     plot( new_horz_line(1,:),new_horz_line(2,:),'r' ,'LineWidth', 1);
%     plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    %set( axis_handle,'NextPlot',hold_state );


end


% DOUG'S ADDITIONS TO FIT_ELLIPSE:
phi = orientation_rad;

% 1. Calculate alpha as the angle between x-axis and the major axis:
if a>b          % a is the major axis, and alpha is the opposite sign as phi
    alpha = -1*phi;
    alpha = [alpha-pi, alpha, alpha+pi];
elseif a<b      % b is the major axis, and alpha = (-1*phi +/- pi/2);
    olda = a;
    a = b;
    b = olda;
    alpha1 = (-1*phi + pi/2);   % could also be alpha = (phi - pi/2).  Depends upon previous values of alpha
    alpha2 = (-1*phi - pi/2);
    alpha = [alpha1-pi, alpha1, alpha1+pi, alpha2-pi, alpha2, alpha2+pi];
end

[minval, minloc] = min(abs(prev_alpha - alpha));

alpha = alpha(minloc);
major = long_axis;
minor = short_axis;
xbar_e = X0_in;
ybar_e = Y0_in;

% AARON'S ADDITIONS TO FIT_ELLIPSE
X_ellipse = rotated_ellipse(1,:);
Y_ellipse = rotated_ellipse(2,:);

% axis equal

end


function [tot_alpha,major,minor] = getrotation(ext_X, ext_Y, prev_alpha)
% Get the accurate angle of rotation for the current cross-section. Used on
% the node and the first 10 cross-sections below the node to verify the
% notch orientation of the stalk

npoints_slice_ext = length(ext_X);

% Uses a fit ellipse function to identify the angle of rotation along the long axis of the cross-section
% (only takes into account the exterior boundaries)
[alpha, ~, ~, ~, ~, ~, ~] = fit_ellipse_R2( ext_X, ext_Y, prev_alpha, gca );

% Reorders and rotates the stalk's exterior and interior
% Rotates an extra 90 degrees so the long axis is vertical
[~, ~, ~, ~, ~, ~, ext_xi, ext_yi, ~, ~] = reorder_V2(ext_X, ext_Y, alpha-pi/2);
% [~, ~, ~, ~, ~, ~, int_xi, int_yi, ~, ~] = reorder_V2_interior(int_X, int_Y, alpha-pi/2, mean(ext_X), mean(ext_Y));

%         close(gcf)
ext_xi = ext_xi';
ext_yi = ext_yi';
% int_xi = int_xi';
% int_yi = int_yi';

% NEW POLAR COORDINATES
ti_ext = 0:2*pi/npoints_slice_ext:2*pi;                             % Creates a theta vector according to the inputted resolution
ti_ext = ti_ext(1:end-1);                                           % The last point is not necessary

for j = 1:length(ti_ext)                                        
    ext_rhoi(:,j) = sqrt(ext_xi(:,j)^2 + ext_yi(:,j)^2);  % Creates an exterior rho vector using the Pythagorean theorem 
end

% LOCATES THE NOTCH -----------------------------------

% Creates the two cut-outs to look for the notch in
window1 = find(ti_ext >   pi/4 & ti_ext < 3*pi/4);
window2 = find(ti_ext > 5*pi/4 & ti_ext < 7*pi/4);

% The 2 windows on all the coordinates
w1_ti = ti_ext(window1);
w2_ti = ti_ext(window2);
w1_rhoi = ext_rhoi(window1);
w2_rhoi = ext_rhoi(window2);
w1_ext_xi = ext_xi(window1);
w2_ext_xi = ext_xi(window2);
w1_ext_yi = ext_yi(window1);
w2_ext_yi = ext_yi(window2);

% Extreme smoothing needed to find peaks
w1_rhoi = smooth(w1_rhoi, 30);
w2_rhoi = smooth(w2_rhoi, 30);

% Peakfinding on the two windows
sel1 = (max(w1_rhoi)-min(w1_rhoi))/32;
w1_peaklocs = peakfinder(w1_rhoi,sel1);
sel2 = (max(w2_rhoi)-min(w2_rhoi))/32;
w2_peaklocs = peakfinder(w2_rhoi,sel2);

% The amount of peaks in each window
w1peaks = length(w1_peaklocs);
w2peaks = length(w2_peaklocs);

if w1peaks >= w2peaks % notch on top
    cut1 = window1(1);
    cut2 = window1(end);
    spin = -pi/2;
elseif w1peaks < w2peaks % notch on bottom
    cut1 = window2(1);
    cut2 = window2(end);
    spin = pi/2;
end

% "Pie" vectors (the external cross sections with the notch cut out)
pier = [ext_rhoi(cut2+1:end)   ext_rhoi(1:cut1-1)];
piet = [ti_ext(cut2+1:end)     ti_ext(1:cut1-1)];
piex = [ext_xi(cut2+1:end)     ext_xi(1:cut1-1)];
piey = [ext_yi(cut2+1:end)     ext_yi(1:cut1-1)];

% Rotate the cross-section again to be horizontal / notch on the right
[~, ~, ~, ~, ~, ~, ext_xi, ext_yi, ~, ~] = reorder_V2(ext_xi, ext_yi, spin);
% [~, ~, ~, ~, ~, ~, int_xi, int_yi, ~, ~] = reorder_V2_interior(int_xi, int_yi, spin, mean(ext_xi), mean(ext_yi));
[~, ~, ~, ~, ~, ~, piex,   piey,   ~, ~] = reorder_V2(piex,   piey,   spin);

% Fitting an ellipse to the cross-section with the notch removed to
% get a more accurate alpha
[new_alpha, major, minor, ~, ~, ~, ~] = fit_ellipse_R2( piex, piey, alpha, gca );
tot_alpha = alpha + new_alpha;

end




function [X, Y, x, y, r, t, xp, yp, rp, tp] = reorder_V2(x, y, alpha)
% function centers and reorders the parametric curves described
% by x and y.
%
% INPUTS:   x,y - column vectors of x & y values
%           alpha - angular orientation of the stalk           
%
% OUTPUTS:  Nomenclature: CAPS variables in image coordinates, lower case relative to stalk coordinate system. "p" indicates "prime" - for a rotated coordinate system.
%           X and Y - reordered x and y values.
%           x, y - reordered x and y values relative to xbar and ybar (but not rotated)
%           r, t - reordered polar coordinates: r and theta values (NOT rotated)
%           xp, yp - reordered x and y values relative to xbar and ybar AND rotated.
%           rp, tp - reordered polar coordinates: r and theta values, AND rotated
%
% Note: designed for use with external contour only (hard to do both since int and ext have different numbers of data points).
%
% VERSION HISTORY
%   V2 - 4/5/18 - updated inputs and comments to reflect changes since conversion from smorder_R2
%
%

% 1. shift origin to center
X = x;
Y = y;

% Locate Center
xbar = mean(x);
ybar = mean(y);

% shift origin to center
y = y - ybar;
x = x - xbar;

% 2. compute theta values in the (xbar, ybar) coordinate system
t = atan2(y,x);        

% 3. correct progression direction    
% make sure points progress in the positive direction - this is done by checking
% the sign of the median difference between theta values. If positive, rotation is positive. If
% negative, all indices are reversed.
theta_diff = diff(t);                   % differences between subsequent theta values
med_diff = median(theta_diff);              % find the MEDIAN theta difference (required because a theta shift value of +/-2pi is always present)

if med_diff < 0                             % check if theta differences are negative
    t = t(end:-1:1);                % if so, reverse indices in all variables
    x = x(end:-1:1);
    y = y(end:-1:1);
    X = X(end:-1:1);
    Y = Y(end:-1:1);
end

t = t(:);
x = x(:);
y = y(:);
X = X(:);
Y = Y(:);

% 4. Reorder variables starting from the value closest to alpha and shift theta values to the stalk coordinate system
temp_theta = t;                         % temp_theta variable used for convenience
index = (t<=alpha);                     % indices of theta values less than alpha
temp_theta(index) = NaN;                    % these values are ignored using NaN
[minval, alphloc] = min(temp_theta);         % find the smallest remaining value
t = [t(alphloc:end); t(1:alphloc-1)];      % all indices are reordered accordingly.
x = [x(alphloc:end); x(1:alphloc-1)];
y = [y(alphloc:end); y(1:alphloc-1)];
X = [X(alphloc:end); X(1:alphloc-1)];
Y = [Y(alphloc:end); Y(1:alphloc-1)];



% 5. Shift theta values so that they are continuous (no jumps) 
% this loop is fairly wasteful, costs about 0.04 seconds.  There are much faster ways, for example just searching for sudden jumps in theta.  But these are likely less robust.
for i = 2:length(t)    
    poss_theta = [t(i) - 2*pi, t(i)-pi, t(i), t(i) + pi, t(i)+2*pi];
    [minval, alphloc] = min(abs(t(i-1) - poss_theta));
    t(i) = poss_theta(alphloc);   
end
% IS THIS STEP EVEN NECESSARY????


% 6. Define radius values
rsq = x.^2 + y.^2;


% 7. Data for coordinate system centered at (xbar, ybar) (no tilt) - First point is at the major axis
x = x; 
y = y;
r = rsq.^0.5;
t = t;

% 8. Data for tilted coordinate system centered at (xbar, ybar) (no tilt)
xp = x*cos(-alpha)-y*sin(-alpha);
yp = x*sin(-alpha)+y*cos(-alpha);
rp = (x.^2 + y.^2).^0.5;
tp = t-alpha;


end




function varargout = peakfinder(x0, sel, thresh, extrema)
%PEAKFINDER Noise tolerant fast peak finding algorithm
%   INPUTS:
%       x0 - A real vector from the maxima will be found (required)
%       sel - The amount above surrounding data for a peak to be
%           identified (default = (max(x0)-min(x0))/4). Larger values mean
%           the algorithm is more selective in finding peaks.
%       thresh - A threshold value which peaks must be larger than to be
%           maxima or smaller than to be minima.
%       extrema - 1 if maxima are desired, -1 if minima are desired
%           (default = maxima, 1)
%   OUTPUTS:
%       peakLoc - The indicies of the identified peaks in x0
%       peakMag - The magnitude of the identified peaks
%
%   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
%       are at least 1/4 the range of the data above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
%       that are at least sel above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local 
%       maxima that are at least sel above surrounding data and larger
%       (smaller) than thresh if you are finding maxima (minima).
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
%       data if extrema > 0 and the minima of the data if extrema < 0
%
%   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
%       local maxima as well as the magnitudes of those maxima
%
%   If called with no output the identified maxima will be plotted along
%       with the input data.
%
%   Note: If repeated values are found the first is identified as the peak
%
% Ex:
% t = 0:.0001:10;
% x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
% x(1250:1255) = max(x);
% peakfinder(x)
%
% Copyright Nathanael C. Yoder 2011 (nyoder@gmail.com)

% Perform error checking and set defaults if not passed in
error(nargchk(1,4,nargin,'struct'));
error(nargoutchk(0,2,nargout,'struct'));

s = size(x0);
flipData =  s(1) < s(2);
len0 = numel(x0);
if len0 ~= s(1) && len0 ~= s(2)
    error('PEAKFINDER:Input','The input data must be a vector')
elseif isempty(x0)
    varargout = {[],[]};
    return;
end
if ~isreal(x0)
    warning('PEAKFINDER:NotReal','Absolute value of data will be used')
    x0 = abs(x0);
end

if nargin < 2 || isempty(sel)
    sel = (max(x0)-min(x0))/4;
elseif ~isnumeric(sel) || ~isreal(sel)
    sel = (max(x0)-min(x0))/4;
    warning('PEAKFINDER:InvalidSel',...
        'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
elseif numel(sel) > 1
    warning('PEAKFINDER:InvalidSel',...
        'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
    sel = sel(1);
end

if nargin < 3 || isempty(thresh)
    thresh = [];
elseif ~isnumeric(thresh) || ~isreal(thresh)
    thresh = [];
    warning('PEAKFINDER:InvalidThreshold',...
        'The threshold must be a real scalar. No threshold will be used.')
elseif numel(thresh) > 1
    thresh = thresh(1);
    warning('PEAKFINDER:InvalidThreshold',...
        'The threshold must be a scalar.  The first threshold value in the vector will be used.')
end

if nargin < 4 || isempty(extrema)
    extrema = 1;
else
    extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
    if extrema == 0
        error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
    end
end

x0 = extrema*x0(:); % Make it so we are finding maxima regardless
thresh = thresh*extrema; % Adjust threshold according to extrema.
dx0 = diff(x0); % Find derivative
dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign

% Include endpoints in potential peaks and valleys
x = [x0(1);x0(ind);x0(end)];
ind = [1;ind;len0];

% x only has the peaks, valleys, and endpoints
len = numel(x);
minMag = min(x);


if len > 2 % Function with peaks and valleys
    
    % Set initial parameters for loop
    tempMag = minMag;
    foundPeak = false;
    leftMin = minMag;
    
    % Deal with first point a little differently since tacked it on
    % Calculate the sign of the derivative since we taked the first point
    %  on it does not neccessarily alternate like the rest.
    signDx = sign(diff(x(1:3)));
    if signDx(1) <= 0 % The first point is larger or equal to the second
        ii = 0;
        if signDx(1) == signDx(2) % Want alternating signs
            x(2) = [];
            ind(2) = [];
            len = len-1;
        end
    else % First point is smaller than the second
        ii = 1;
        if signDx(1) == signDx(2) % Want alternating signs
            x(1) = [];
            ind(1) = [];
            len = len-1;
        end
    end
    
    % Preallocate max number of maxima
    maxPeaks = ceil(len/2);
    peakLoc = zeros(maxPeaks,1);
    peakMag = zeros(maxPeaks,1);
    cInd = 1;
    % Loop through extrema which should be peaks and then valleys
    while ii < len
        ii = ii+1; % This is a peak
        % Reset peak finding if we had a peak and the next peak is bigger
        %   than the last or the left min was small enough to reset.
        if foundPeak
            tempMag = minMag;
            foundPeak = false;
        end
        
        % Make sure we don't iterate past the length of our vector
        if ii == len
            break; % We assign the last point differently out of the loop
        end
        
        % Found new peak that was lager than temp mag and selectivity larger
        %   than the minimum to its left.
        if x(ii) > tempMag && x(ii) > leftMin + sel
            tempLoc = ii;
            tempMag = x(ii);
        end
        
        ii = ii+1; % Move onto the valley
        % Come down at least sel from peak
        if ~foundPeak && tempMag > sel + x(ii)
            foundPeak = true; % We have found a peak
            leftMin = x(ii);
            peakLoc(cInd) = tempLoc; % Add peak to index
            peakMag(cInd) = tempMag;
            cInd = cInd+1;
        elseif x(ii) < leftMin % New left minima
            leftMin = x(ii);
        end
    end
    
    % Check end point
    if x(end) > tempMag && x(end) > leftMin + sel
        peakLoc(cInd) = len;
        peakMag(cInd) = x(end);
        cInd = cInd + 1;
    elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
        peakLoc(cInd) = tempLoc;
        peakMag(cInd) = tempMag;
        cInd = cInd + 1;
    end
    
    % Create output
    peakInds = ind(peakLoc(1:cInd-1));
    peakMags = peakMag(1:cInd-1);
else % This is a monotone function where an endpoint is the only peak
    [peakMags,xInd] = max(x);
    if peakMags > minMag + sel
        peakInds = ind(xInd);
    else
        peakMags = [];
        peakInds = [];
    end
end

% Apply threshold value.  Since always finding maxima it will always be
%   larger than the thresh.
if ~isempty(thresh)
    m = peakMags>thresh;
    peakInds = peakInds(m);
    peakMags = peakMags(m);
end



% Rotate data if needed
if flipData
    peakMags = peakMags.';
    peakInds = peakInds.';
end



% Change sign of data if was finding minima
if extrema < 0
    peakMags = -peakMags;
    x0 = -x0;
end
% Plot if no output desired
if nargout == 0
    if isempty(peakInds)
        disp('No significant peaks found')
    else
        figure;
        plot(1:len0,x0,'b',peakInds,peakMags,'ro','linewidth',2);
    end
else
    varargout = {peakInds,peakMags};
end

end