function fit_and_refit_ellipse(FileName,t,SaveName)
% Fit an ellipse to a given cross section, then find the notch, cut it out
% of the data, and refit the ellipse. Finally, offset the ellipse with a
% constant rind thickness.
close;

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'ext_T','ext_R','int_T','int_R');

% Convert from polar to Cartesian coordinates for the ellipse fit function
[ext_x,ext_y] = pol2cart(ext_T,ext_R);
% [int_x,int_y] = pol2cart(int_T,int_R);

npoints = 360;

%% EXTERIOR BOUNDARIES
npoints_slice = length(ext_x);

% Uses a fit ellipse function to identify the angle of roatation along the long axis of the cross-section
alpha = 0;
prev_alpha = 0;
[alpha, ellx, elly, major, minor, ~, ~] = fit_ellipse_R3(ext_x, ext_y, prev_alpha, gca);
alpha
% major
% minor
[Tnew,Rnew] = cart2pol(ellx,elly);

figure(1)
polarplot(ext_T,ext_R);
hold on
polarplot(Tnew,Rnew);


% Reorders and rotates the stalk's exterior
% Rotates an extra 90 degrees so the long axis is veritcal
[~, ~, ~, ~, ~, ~, ext_xi, ext_yi, ~, ~] = reorder_V2(ext_x, ext_y, pi/2);
% figure(2)
% plot(ext_xi,ext_yi)

% close(gcf)
ext_xi = ext_xi';
ext_yi = ext_yi';
% 
% NEW POLAR COORDINATES
ti = 0:2*pi/npoints_slice:2*pi;                             % Creates a theta vector according to the inputted resolution
ti = ti(1:end-1);                                           % The last point is not necessary

for j = 1:length(ti)                                        
    rhoi(:,j) = sqrt(ext_xi(:,j)^2 + ext_yi(:,j)^2);  % Creates a rho vector using pythagorean's theorem 
end
% 
% % LOCATES THE NOTCH -----------------------------------
% 
% Creates the two cut-outs to look for the notch in
window1 = find(ti >   pi/4 & ti < 3*pi/4);
window2 = find(ti > 5*pi/4 & ti < 7*pi/4);

% The 2 windows on all the coordinates
w1_ti = ti(window1);
w2_ti = ti(window2);
w1_rhoi = rhoi(window1);
w2_rhoi = rhoi(window2);
% w1_ext_xi = ext_xi(window1);
% w2_ext_xi = ext_xi(window2);
% w1_ext_yi = ext_yi(window1);
% w2_ext_yi = ext_yi(window2);

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

% Force the variables into rows
rhoi = rhoi(:)';
ti = ti(:)';
% ext_xi = ext_xi(:)';
% ext_yi = ext_yi(:)';

% "Pie" vectors (the cross sections with the notch cut out)
pier = [rhoi(cut2+1:end)   rhoi(1:cut1-1)];
piet = [ti(cut2+1:end)     ti(1:cut1-1)];
% figure(3)
% polarplot(piet,pier)
% piex = [ext_xi(cut2+1:end) ext_xi(1:cut1-1)];
% piey = [ext_yi(cut2+1:end) ext_yi(1:cut1-1)];
% 
% 
% % Rotate the cross-section again to be horizontal / notch on the right
% [~, ~, ~, ~, ~, ~, ext_xi, ext_yi, ~, ~] = reorder_V2(ext_xi, ext_yi, spin);
% [~, ~, ~, ~, ~, ~, piex,   piey,   ~, ~] = reorder_V2(piex,   piey,   spin);
% 


[piex,piey] = pol2cart(piet,pier);
% Fitting an ellipse to the cross-section with the notch removed to
% get a more accurate alpha
[alpha, ellx_refit, elly_refit, major_refit, minor_refit, ~, ~] = fit_ellipse_R3(piex, piey, alpha, gca);
alpha
% major_refit
% minor_refit
[Tnew,Rnew] = cart2pol(ellx_refit,elly_refit);

figure(1);
polarplot(Tnew,Rnew);
% 
% % Rotating according to the new, more accurate alpha
% [~, ~, ~, ~, ~, ~, ext_xi, ext_yi, ~, ~] = reorder_V2(ext_xi, ext_yi, new_alpha); 
% 
% % % Scaling by the major and minor axes
% % ext_xscaled = ext_xi ./ ((major + minor)/2);
% % ext_yscaled = ext_yi ./ ((major + minor)/2);
% % 
% % % Downsampling
% % idx =  1:length(ext_xi);                               % Index
% % idxq = linspace(min(idx), max(idx), npoints);               % Interpolation Vector
% % ext_xi = interp1(idx, ext_xi, idxq, 'linear');         % Downsampled Vector
% % 
% % idx = 1:length(ext_yi);                                % Index
% % idxq = linspace(min(idx), max(idx), npoints);               % Interpolation Vector
% % ext_yi = interp1(idx, ext_yi, idxq, 'linear');         % Downsampled Vector
% 
% 
% % Getting all the first indicies to be exactly on the x-axis
% m_ext = (ext_yi(1)-ext_yi(end))/(ext_xi(1)-ext_xi(end));    % Solve for slope
% ext_x1 = ext_xi(1);                                         % x-point on the right line
% ext_y1 = ext_yi(1);                                         % y-point on the right line
% b_ext = ext_y1 - m_ext*ext_x1;                              % Solve for y-intercept
% ext_y2 = 0;                                                 % We want the x-value where y=0
% ext_x2 = (ext_y2-b_ext)/m_ext;                              % Solving for the x-value
% 
% % Shift all the data according the the differences
% xdif1 = ext_xi(1) - ext_x2;
% ext_xi = ext_xi - xdif1; 
% ydif1 = ext_yi(1);
% ext_yi = ext_yi - ydif1; 
% 
% % NEW NEW POLAR COORDINATES
% tDCSR = 0:2*pi/npoints:2*pi;                             % Creates a theta vector according to the inputted resolution
% tDCSR = tDCSR(1:end-1);                                           % The last point is not necessary
% 
% for j = 1:length(tDCSR)                                        
%     rhoDCSR(j) = sqrt(ext_xi(j)^2 + ext_yi(j)^2);  % Creates a rho vector using pythagorean's theorem 
% end
% 
% polarplot(tDCSR,rhoDCSR)

% Save new ext_T and ext_R for output mat file
ext_T = Tnew';
ext_R = Rnew';
int_T = Tnew';
int_R = zeros(size(ext_R));

% Offset ext_R by t to get a very rough elliptical fit
for i = 1:length(ext_R)
    int_R(i) = ext_R(i) - t;    
end


% Save the final data in a new mat file
SaveFile       = fullfile(FolderName, SaveName);
save(SaveFile,'ext_T','ext_R','int_T','int_R');
end