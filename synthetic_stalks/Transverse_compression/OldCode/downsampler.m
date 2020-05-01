function [r_sample, theta_sample, x_sample, y_sample] = downsampler(r, t, xbar, ybar, N)

% function to perform downsampling
%
%   INPUTS: rf, thetaf - final (i.e. smoothed) r and theta data sets
%           xbar, ybar - center points of x and y curves
%           N - number of points to extract between theta = 0 and theta =
%           360 degrees. Point #1 is at theta = 0.
%
%   OUTPUTS:    sampled r, theta, x, and y data.

theta_sample = [0:2*pi/N:2*pi-2*pi/N]';                              % N points, starting at theta = 0;
r_sample = interp1(t,r,theta_sample,'linear', 'extrap');      % cubic interpolation, extrapolation needed for theta = 0

x_sample = r_sample.*cos(theta_sample) + xbar;                      % converting to x and y coordinate system
y_sample = r_sample.*sin(theta_sample) + ybar;

% % Visualization - plotting original data, smoothed curves, and downsampled data
% plot(x,y,'k')
% hold on;
% plot(xf,yf,'b', 'LineWidth', 2)
% plot(x_sample,y_sample, '.r')

