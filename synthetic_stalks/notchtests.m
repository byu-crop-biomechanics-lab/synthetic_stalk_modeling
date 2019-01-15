clear all;
close all;

dmaj = 2;
dmin = 1;

nwidth = 0.5;
ndepth = 0.1;
if nwidth >= dmaj
    error('Notch width is too large. Increase dmaj or decrease nwidth.')
end

theta = linspace(-pi,pi,500);   % The range of theta must stay as -pi to pi
yshift = zeros(1,length(theta));
notch = zeros(1,length(theta));
y = (dmin/2)*cos(theta);
x = (dmaj/2)*sin(theta);
% flatlimit = pi/4;
flatlimit = asin(nwidth/dmaj);
for i = 1:length(notch)
    if theta(i) < -flatlimit | theta(i) > flatlimit
        notch(i) = 0;
    else
        notch(i) = ((ndepth/2)*(cos((2*pi/nwidth)*theta(i)))+(ndepth/4));
    end
end

% flatheight = 0;

for i = 1:length(yshift)
    % Don't shift y if outside the chosen bounds
    if theta(i) < -pi | theta(i) > pi
        yshift(i) = 0;
    elseif theta(i) >= -flatlimit & theta(i) <= flatlimit
%         yshift(i) = (dmin/2)*cos(-flatlimit);   % Keep the same value as the value at flatlimits
        yshift(i) = y(i) - notch(i);
    else
        yshift(i) = y(i);
    end          
end

% yshift(yshift == 0) = NaN;  % THIS STILL CAUSES A HOLE IN THE BOTTOM OF THE ELLIPSE
% for i = 1:length(theta)
%     if theta(i) == 0
%         yshift(i) = 0;
%     end
% end
% x(1) = 0;
% x(end) = 0;

% figure(1)
% plot(theta,y)
% hold on
% plot(theta,yshift)
% legend('y','yshift');
% xlabel('\theta, radians')

hold off
figure(2)
plot(x,yshift)

figure(3)
plot(theta,notch)