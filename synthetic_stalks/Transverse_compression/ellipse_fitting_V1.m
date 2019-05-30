function ellipse_fitting_V1(FileName,SaveName)
% Load a mat file that has exterior XY data and avgrindthickness data.
% Cycle through cross sections and select the angular range that
% contains the notch so it's ignored during ellipse fitting. Then save the
% ellipses in a mat file. Take the difference between the interior and
% exterior points and their ellipse approximations and save those as well
% for later PCA.

close all

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'avg_rind_thick','ext_rhoDCSR','ext_tDCSR','ext_xDCSR','ext_yDCSR',...
    'int_rhoDCSR','int_tDCSR','int_xDCSR','int_yDCSR');

% Make copies of original data to work with
R_ext = ext_rhoDCSR;
R_int = int_rhoDCSR;
T = (ext_tDCSR(:,1));
X_ext = ext_xDCSR;
Y_ext = ext_yDCSR;
X_int = int_xDCSR;
Y_int = int_yDCSR;

% % Plot in Cartesian coordinates to check results
% plot(X_ext(:,1),Y_ext(:,1),'.','LineWidth',2);
% hold on
% plot(X_int(:,1),Y_int(:,1),'.','LineWidth',2);
% pause();
% close;


N = length(avg_rind_thick);

A = zeros(N,1);
B = zeros(N,1);
ALPHA = zeros(N,1);


ELLIPSE_XY = zeros(N,360,2);
ELLIPSE_CENTERS = zeros(N,2);
ELLIPSE_T = zeros(N,360);
ELLIPSE_R_ext = zeros(N,360);
ELLIPSE_R_int = zeros(N,360);

DIFF_R_ext = zeros(N,360);
DIFF_R_int = zeros(N,360);
R_ext = zeros(N,360);
R_int = zeros(N,360);

AVG_RIND_T = zeros(N,1);

prev_alpha = 0;

min_angle = 135;
min_angle = min_angle*(pi/180);     % Convert angle to radians
max_angle = 225;
max_angle = max_angle*(pi/180);

for i = 1:N
    % Define the notch range
    for j = 1:360
        if T(j) > min_angle
            min_index = j-1;
            break
        end
    end
    
    for j = 1:360
        if T(j) > max_angle
            max_index = j;
            break
        end
    end
    
    % Cut out the notch from the XY data
    window = [linspace(1,min_index,min_index),linspace(max_index,360,(360-max_index + 1))];
    xcut = X_ext(window,i);
    ycut = Y_ext(window,i);
%     figure(2);
%     plot(X_ext(:,i),Y_ext(:,i));
%     hold on
%     plot(xcut,ycut);
%     pause();
%     close;

    [alpha, major, minor, xbar_e, ybar_e, X_ellipse, Y_ellipse] = fit_ellipse_R2( xcut, ycut, prev_alpha, gca );
%     alpha
%     xbar_e
%     ybar_e
%     major
%     minor
    
    
    % NEED TO SHIFT ELLIPSE CENTER TO THE GEOMETRIC CENTER IN ORDER TO
    % COMPUTE THE TRANSFORMATION TO POLAR COORDINATES AS BELOW!
    
    % Save ellipse center shift
    ELLIPSE_CENTERS(i,:) = [mean(X_ellipse), mean(Y_ellipse)];
    
    % Shift XY data before converting to polar
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
    
%     % Plot in Cartesian coordinates to check results
%     plot(X_ext_shift,Y_ext_shift,'.','LineWidth',2);
%     hold on
%     plot(X_int_shift,Y_int_shift,'.','LineWidth',2);
%     pause();
%     close;
    
    
    % Convert X_ellipse and Y_ellipse to polar coordinates
    theta = 0:2*pi/360:2*pi;
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
    ELLIPSE_R_int(i,:) = ext_rho_ellipse - avg_rind_thick(i);
    R_ext(i,:) = ext_rho;
    R_int(i,:) = int_rho;
    
    % Get difference between the ellipse and the real data (if the ellipse
    % overestimates, then the value of DIFF will be positive)
    DIFF_R_ext(i,:) = ext_rho_ellipse - ext_rho;
    DIFF_R_int(i,:) = int_rho_ellipse - int_rho;
    
%     DIFF_R_ext(i,:) = ext_rho_ellipse - transpose(R_ext(:,i));
%     DIFF_R_int(i,:) = int_rho_ellipse - transpose(R_int(:,i));
    
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
