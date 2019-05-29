function ellipse_fitting(FileName,N,SaveName)
% Load a mat file that has exterior XY data and avgrindthickness data.
% Choose five random cross sections and select the angular range that
% contains the notch so it's ignored during ellipse fitting. Then save the
% ellipses in a mat file. Take the difference between the interior and
% exterior points and their ellipse approximations and save those as well
% for later PCA.

close all

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'ext_rhoDCSR','ext_xDCSR','ext_yDCSR','int_rhoDCSR','int_xDCSR',...
    'int_yDCSR','tDCSR','avg_rind_thickness');

indices = randi([1 360],N,1);
if length(indices) ~= length(unique(indices))
    error('An index was repeated');
end

A = zeros(N,1);
B = zeros(N,1);
ALPHA = zeros(N,1);

ELLIPSE_XY = zeros(N,360,2);
ELLIPSE_T = zeros(N,360);
ELLIPSE_R = zeros(N,360);

DIFF_ext_R = zeros(N,360);
DIFF_int_R = zeros(N,360);

AVG_RIND_T = zeros(N,1);

prev_alpha = 0;
prompt1 = 'Enter the minimum angle that defines the notch, in degrees.';
prompt2 = 'Enter the maximum angle that defines the notch, in degrees.';

for i = 1:N
    indices(i)
    polarplot(tDCSR,ext_rhoDCSR(1,:,indices(i)));
    hold on
    polarplot(tDCSR,int_rhoDCSR(1,:,indices(i)));
    
    % Visually determine the angular range that defines the notch
    min_angle = input(prompt1);
    min_angle = min_angle*(pi/180);     % Convert angle to radians
    pause();
    max_angle = input(prompt2);
    max_angle = max_angle*(pi/180);
    close;
    
    % Define the notch range
    for j = 1:360
        if tDCSR(j) > min_angle
            min_index = j-1;
            break
        end
    end
    
    for j = 1:360
        if tDCSR(j) > max_angle
            max_index = j;
            break
        end
    end
    
    % Cut out the notch from the XY data
    window = [linspace(1,min_index,min_index),linspace(max_index,360,(360-max_index + 1))];
    xcut = ext_xDCSR(1,window,indices(i));
    ycut = ext_yDCSR(1,window,indices(i));
    figure(2);
    plot(ext_xDCSR(1,:,indices(i)),ext_yDCSR(1,:,indices(i)));
    hold on
    plot(xcut,ycut);
    pause();
    close;
    [alpha, ellx, elly, major, minor, ~, ~] = fit_ellipse_R3(xcut, ycut, prev_alpha, 360, gca);
    
    % Convert ellx and elly to polar coordinates
    theta = 0:2*pi/360:2*pi;
    theta = theta(1:end-1);
    ext_rho = zeros(size(theta));
    for j = 1:length(theta)
        ext_rho(j) = sqrt(ellx(j)^2 + elly(j)^2);
    end
    
    % Interpolate ellipse points so they are at the same theta values as
    % tDCSR
    rho_new = interp1(theta,ext_rho,tDCSR);
    assignin('base','theta',theta);
    ext_rho = rho_new;
    assignin('base','ext_rho',ext_rho);
    int_rho = ext_rho - avg_rind_thickness(indices(i));
    theta = tDCSR;
    
    polarplot(tDCSR,ext_rhoDCSR(1,:,indices(i)));
    hold on
    polarplot(tDCSR,int_rhoDCSR(1,:,indices(i)));
    polarplot(theta,ext_rho);
    polarplot(theta,int_rho);
    pause();
    close;
    
    A(i) = major;
    B(i) = minor;
    ALPHA(i) = alpha;

    ELLIPSE_XY(i,:,1) = ellx;
    ELLIPSE_XY(i,:,2) = elly;
    ELLIPSE_T(i,:) = theta;
    ELLIPSE_R(i,:) = ext_rho;
    
    % Get difference between the ellipse and the real data (if the ellipse
    % overestimates, then the value of DIFF will be positive)
    DIFF_ext_R(i,:) = ext_rho - ext_rhoDCSR(1,:,indices(i));
    DIFF_int_R(i,:) = int_rho - int_rhoDCSR(1,:,indices(i));
    
    AVG_RIND_T(i) = avg_rind_thickness(indices(i));
    
end

% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R','DIFF_ext_R','DIFF_int_R','AVG_RIND_T','indices');

end


