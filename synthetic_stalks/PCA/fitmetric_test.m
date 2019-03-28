% fitmetric_test.m: Testing of fit for synthetic points against the
% original data, using a circle.

close;
clear;

N = 100;
theta = linspace(0,2*pi,N);
circlex = cos(theta);
circley = sin(theta);

% Noisy data approximation case
noisex = unifrnd(-0.05,0.05,1,N);
noisey = unifrnd(-0.05,0.05,1,N);
fitnoisex = circlex + noisex;
fitnoisey = circley + noisey;

% All values are off by the same radial amount case
bigx = 1.25*circlex;
bigy = 1.25*circley;

plot(circlex,circley)
daspect([1 1 1])
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);

hold on
plot(bigx,bigy)
plot(fitnoisex,fitnoisey)

legend('Original data: perfect circle','Noisy fit','Perfect circle off by set amount');

noisefit = getfit(circlex,fitnoisex,circley,fitnoisey)
bigfit = getfit(circlex,bigx,circley,bigy)

% function [fitmetric] = getfit(xreal,yreal,xsynth,ysynth)
%     % Get the overall fit of a synthetic shape against a real shape, in a
%     % single value similar to a standard deviation
%     Nxreal = length(xreal);
%     Nyreal = length(yreal);
%     Nxsynth = length(xsynth);
%     Nysynth = length(ysynth);
%     
%     assert(Nxreal == Nxsynth);
%     assert(Nyreal == Nysynth);
%     
%     sum = 0;
%     
%     for i = 1:Nxreal
%         xsq = (xsynth(i) - xreal(i))^2;
%         ysq = (ysynth(i) - yreal(i))^2;
%         sum = sum + xsq + ysq;
%     end
%     
%     fitmetric = sqrt(sum/Nxreal);
%     
% end