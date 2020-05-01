function ellipse_fitting_opt(FileName,SaveName)
% Load a mat file that has exterior XY data and avgrindthickness data.
% Cycle through cross sections and select the angular range that
% contains the notch so it's ignored during ellipse fitting. Then save the
% ellipses in a mat file. Take the difference between the interior and
% exterior points and their ellipse approximations and save those as well
% for later PCA.

close all

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'avg_rind_thickness','R_ext','R_int','T','X_ext','X_int',...
    'Y_ext','Y_int');

N = length(avg_rind_thickness);

A = zeros(N,1);
B = zeros(N,1);
% ALPHA = zeros(N,1);

ELLIPSE_XY = zeros(N,360,2);
ELLIPSE_T = zeros(N,360);
ELLIPSE_R_ext = zeros(N,360);
ELLIPSE_R_int = zeros(N,360);

DIFF_R_ext = zeros(N,360);
DIFF_R_int = zeros(N,360);

AVG_RIND_T = zeros(N,1);

prev_alpha = 0;
% prompt1 = 'Enter the minimum angle that defines the notch, in degrees.';
% prompt2 = 'Enter the maximum angle that defines the notch, in degrees.';

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
    Rcut = R_ext(window,i);
    Tcut = T(window);
%     figure(2);
%     polarplot(T,R_ext(:,i));
%     hold on
%     polarplot(Tcut,Rcut);
%     pause();
%     close;

    [xopt, ~, ~, ~] = fit_ellipse_opt(Tcut,Rcut);
    major = xopt(1);
    minor = xopt(2);
    
    ext_rho = rpts_ellipse(T,major,minor);
    major_int = major - 2*avg_rind_thickness(i);
    minor_int = minor - 2*avg_rind_thickness(i);
    int_rho = rpts_ellipse(T,major_int,minor_int);
    theta = T;
    
%     % Convert ellx and elly to polar coordinates
%     theta = 0:2*pi/360:2*pi;
%     theta = theta(1:end-1);
%     ext_rho = zeros(size(theta));
%     for j = 1:length(theta)
%         ext_rho(j) = sqrt(ellx(j)^2 + elly(j)^2);
%     end
%     
%     % Interpolate ellipse points so they are at the same theta values as
%     % T
%     rho_new = interp1(theta,ext_rho,T);
%     assignin('base','theta',theta);
%     ext_rho = rho_new;
%     assignin('base','ext_rho',ext_rho);
%     int_rho = ext_rho - avg_rind_thickness(i);
%     theta = T;
    
%     polarplot(T,R_ext(:,i));
%     hold on
%     polarplot(T,R_int(:,i));
%     polarplot(theta,ext_rho,'k','LineWidth',2);
%     polarplot(theta,int_rho,'k','LineWidth',2);
%     pause();
%     close;
    
    A(i) = major;
    B(i) = minor;

    ELLIPSE_T(i,:) = theta;
    ELLIPSE_R_ext(i,:) = ext_rho;
    ELLIPSE_R_int(i,:) = ext_rho - avg_rind_thickness(i);
    
    % Get difference between the ellipse and the real data (if the ellipse
    % overestimates, then the value of DIFF will be positive)
    DIFF_R_ext(i,:) = transpose(ext_rho' - R_ext(:,i));
    DIFF_R_int(i,:) = transpose(int_rho' - R_int(:,i));
    
    AVG_RIND_T(i) = avg_rind_thickness(i);
    
end

% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int',...
    'DIFF_R_ext','DIFF_R_int','AVG_RIND_T','R_ext','R_int','T');

end


function [r] = rpts(N,theta,dmaj,dmin)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end

