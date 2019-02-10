% Generate lots of stalk cross sections with random noise and feature
% variation. 
clear all;
close all;
clc;

% Choose the number of data points to define the stalk shape
N = 360;
theta = linspace(0,2*pi,N);

% Choose how many stalk cross sections to generate:
n = 50;

% Create an empty array (n x N x 2) to represent the x and y data for all
% of the cross sections (a row in slice 1 represents x, and a row in slice
% 2 represents y).
sections = zeros(n,N,2);

% Choose size limits for major and minor diameter (lower bound on major
% diameter must be greater than the upper bound for minor diameter):
dmin_low = 15;
dmin_up = 20;
dmaj_low = dmin_up;
dmaj_up = 25;

aAmplim = 0.125;

%% Main loop
for i = 1:n
    dmaj = unifrnd(dmaj_low,dmaj_up);
    dmin = unifrnd(dmin_low,dmin_up);
    ndepth = unifrnd(0.1,0.75);
    nwidth = unifrnd(1,6);
    nloc = unifrnd(pi-0.2,pi+0.2);
    rotate_angle = unifrnd(-pi,pi);
    
    xasymmetry = Asymmetry(aAmplim,theta,N);
    yasymmetry = Asymmetry(aAmplim,theta,N);
    
    xshift = unifrnd(-10,10);
    yshift = unifrnd(-10,10);
    
    % Random Noise
    noisex = unifrnd(-0.0025,0.0025,1,N);
    noisey = unifrnd(-0.0025,0.0025,1,N);
    
    notch = notch_fn(N,ndepth,nwidth,nloc,theta);
    
    x = xpts(N,theta,notch,dmaj,noisex,xasymmetry);
    y = ypts(N,theta,dmin,noisey,yasymmetry);
    
    [x,y] = rotate(x,y,rotate_angle,N);
    
    x = xshift + x;
    y = yshift + y;
    
    sections(i,:,1) = x;
    sections(i,:,2) = y;
    
end

%% Plot cross sections to verify that none of them are unrealistic
for i = 1:n
    plot(sections(i,:,1),sections(i,:,2));
    hold on
    pause(0.25);    
end


%% Functions
function [x] = xpts(N,theta,notch,dmaj,noisex,asymmetry)
    x = zeros(1,N);
    for i = 1:N
        x(i) = dmaj*(cos(theta(i)) + notch(i) + noisex(i) + asymmetry(i));
    end
end

function [y] = ypts(N,theta,dmin,noisey,asymmetry)
    y = zeros(1,N);
    for i = 1:N
        y(i) = dmin*(sin(theta(i)) + noisey(i) + asymmetry(i));
    end
end

function [notch] = notch_fn(N,ndepth,nwidth,nloc,theta)
    notch = zeros(1,N);
    for i = 1:N
        notch(i) = ndepth/cosh((10/nwidth)*(theta(i)-nloc))^2;
    end
end

function [asymmetry] = Asymmetry(aAmplim,theta,N)
    asymmetry = zeros(1,N);
    aAmp = unifrnd(-aAmplim,aAmplim);
    aSym = unifrnd(-pi,pi);
    asymmetry = aAmp*sin(theta - aSym);
end

function [xrotate,yrotate] = rotate(x,y,psi,N)
    xrotate = zeros(1,N);
    yrotate = zeros(1,N);
    
    R = [cos(psi) -sin(psi); sin(psi) cos(psi)];
    
    for i = 1:N
        temp = [x(i);y(i)];
        temp = R*temp;
        xrotate(i) = temp(1);
        yrotate(i) = temp(2);
    end
end