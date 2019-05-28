%% stalk_cross_sections_polar.m
% RL - 4/11/2019
% Polar implementation of stalk_cross_sections.m
% Generates lots of stalk cross sections with random noise and feature
% variation. 
clear all;
close all;
clc;

% Choose the number of data points to define the stalk shape
N = 180;
theta = linspace(0,2*pi,N);

% Choose how many stalk cross sections to generate:
n = 100;

% Create an empty array (n x N x 2) to represent the x and y data for all
% of the cross sections (a row in slice 1 represents x, and a row in slice
% 2 represents y).
sections = zeros(n,N);

% Choose size limits for major and minor diameter (lower bound on major
% diameter must be greater than the upper bound for minor diameter):
dmin_low = 15;
dmin_up = 20;
dmaj_low = dmin_up;
dmaj_up = 25;

aAmplim = 5;

%% Main loop
for i = 1:n
    dmaj = unifrnd(dmaj_low,dmaj_up);
    dmin = unifrnd(dmin_low,dmin_up);
    
    ndepth = unifrnd(1,5);
    nwidth = unifrnd(1,5);
    nloc = unifrnd(pi-0.2,pi+0.2);      
    
    asymmetry = Asymmetry(aAmplim,theta);    
    
    % Random noise in shape to prevent them from being perfectly smooth
    noise = unifrnd(-0.05,0.05,1,N);
    
    % Generate points to define the cross section
    notch = notch_fn(N,ndepth,nwidth,nloc,theta);    
    r = rpts(N,theta,dmaj,dmin,noise,asymmetry,notch);

%     % Scale the x and y points by a factor related to dmin and dmaj
%     factor = 1/(dmaj + dmin);
% %     factor = 1;
%     x = x*factor;
%     y = y*factor;
    
    % Place cross section data in the larger array of data
    sections(i,:) = r;

    
end

% % Plot cross sections to verify that they're realistic enough
% for i = 1:n
%     polarplot(sections(i,:));
% %     hold on
%     pause(0.1);
% end

% Save data as a mat file for ease of use
save polar_cross_sections.mat sections





%% Functions
function [r] = rpts(N,theta,dmaj,dmin,noise,asymmetry,notch)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2) + noise(i) + asymmetry(i) - notch(i);
    end
end

function [notch] = notch_fn(N,ndepth,nwidth,nloc,theta)
    notch = zeros(1,N);
    for i = 1:N
        notch(i) = ndepth/cosh((10/nwidth)*(theta(i)-nloc))^2;
    end
end

function [asymmetry] = Asymmetry(aAmplim,theta)
%     asymmetry = zeros(1,N);
    aAmp = unifrnd(-aAmplim,aAmplim);
    aSym = unifrnd(-pi,pi);
    asymmetry = aAmp*sin(theta - aSym);
end