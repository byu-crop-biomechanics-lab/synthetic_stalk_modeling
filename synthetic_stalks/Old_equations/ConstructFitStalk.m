%% ConstructFitStalk.m
% Take the fit parameter values from the fit optimization in
% fitcheck_polar.m and construct an STL (requires conversion of points from
% polar to Cartesian in order to work with surf2stl)

% Load parameter values
clear;
load RealStalkFit.mat

nsections = 990;
Npts = 90;
sections = zeros(nsections,Npts);

theta = linspace(0,2*pi,Npts);

% Iterate through the cross sections
for i = 1:nsections
    % Set parameter values for current cross section
    dmaj =  dmajvals(i);
    dmin =  dminvals(i);
    ndepth =  ndepthvals(i);
    nwidth = nwidthvals(i);
    nloc = nlocvals(i);
    aAmp = aAmpvals(i);
    aSym = aSymvals(i);
    
    % Calculate cross section shape
    asymmetry = aAmp*sin(theta - aSym);
    notch = notch_fn(Npts,ndepth,nwidth,nloc,theta);
    r = rpts(Npts,theta,dmaj,dmin,asymmetry,notch);
    
    sections(i,:) = r;
    
end

% % Verify the fit cross sections are being read in correctly
% for i = 1:nsections
%     polarplot(theta,sections(i,:));
%     pause(0.1);
% end

% Convert from polar to Cartesian coordinates
X = zeros(nsections,Npts);
Y = zeros(nsections,Npts);

for i = 1:nsections
    for j = 1:Npts
        X(i,j) = sections(i,j)*cos(theta(j));
        Y(i,j) = sections(i,j)*sin(theta(j));
    end
end

% % Verify the conversion by plotting the cross sections
% for i = 1:nsections
%     plot(X(i,:),Y(i,:));
%     pause(0.1);    
% end

% Create Z data
Z = zeros(nsections,Npts);
a = 0.1;    % Arbitrary distance between cross sections
for i = 1:nsections
    Z(i,:) = i*a;
end

surf2stl_V1('FitStalk.stl',X,Y,Z)

% ------------Stalk geometry functions---------------
function [r] = rpts(N,theta,dmaj,dmin,asymmetry,notch)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2) + asymmetry(i) - notch(i);
    end
end

function [notch] = notch_fn(N,ndepth,nwidth,nloc,theta)
    notch = zeros(1,N);
    for i = 1:N
        notch(i) = ndepth/cosh((10/nwidth)*(theta(i)-nloc))^2;
    end
end