clear all;
close all;

dmaj = 3;
dmin = 2;

nwidth = 2.5;

theta = linspace(-pi,pi,100);
yshift = zeros(1,length(theta));
y = (dmin/2)*cos(theta);
x = (dmaj/2)*sin(theta);
% flatlimit = pi/4;
flatlimit = asin(nwidth/dmaj);
% flatheight = 0;

for i = 1:length(yshift)
    % Don't shift y if outside the chosen bounds
    if theta(i) < -pi | theta(i) > pi
        yshift(i) = 0;
    elseif theta(i) >= -flatlimit & theta(i) <= flatlimit
        yshift(i) = (dmin/2)*cos(-flatlimit);
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

figure(1)
plot(theta,y)
hold on
plot(theta,yshift)

hold off
figure(2)
plot(x,yshift)