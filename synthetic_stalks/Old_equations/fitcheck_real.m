% Visually check the fit produced by stalk_cross_fit.m function
load AAbounds.mat ext_X ext_Y
numsection = 1;
xdata = ext_X;
ydata = ext_Y;

[xopt, fopt, exitflag, output] = stalk_cross_fit_real(xdata,ydata);

% load XYryan.mat ext_xDCSR ext_yDCSR
% 
% numsection = 1; % MAKE SURE THIS IS THE SAME CROSS SECTION THAT WAS OPTIMIZED
% 
% xdata = ext_xDCSR(1,:,numsection);
% ydata = ext_yDCSR(1,:,numsection);

% Get variable values from the output of stalk_cross_fit.m:
dmaj            = xopt(1);
dmin            = xopt(2); 
ndepth          = xopt(3);
nwidth          = xopt(4);
nloc            = xopt(5);
rotate_angle    = xopt(6);
xshift          = xopt(7);
yshift          = xopt(8);
xaAmp           = xopt(9);
xaSym           = xopt(10);
yaAmp           = xopt(11);
yaSym           = xopt(12);

N = length(xdata);
theta = linspace(0,2*pi,N);
phi = nloc - pi;

% Analysis functions
xasymmetry = xaAmp*sin(theta - xaSym);
yasymmetry = yaAmp*sin(theta - yaSym);
xnotch = notch_fn(N,ndepth,nwidth,nloc,theta);
ynotch = zeros(1,N);
[xnotch_rot,ynotch_rot] = rotate(xnotch,ynotch,phi,N);

xmain = xpts(N,theta,dmaj,xasymmetry);
ymain = ypts(N,theta,dmin,yasymmetry);

x = xmain + xnotch_rot;
y = ymain + ynotch_rot;
[x,y] = rotate(x,y,rotate_angle,N);
xsynth = xshift + x;
ysynth = yshift + y;



% Plot the original data and the fit data on top of each other
plot(xdata,ydata);
hold on
plot(xsynth,ysynth)
legend('Original data','Fit model');





% ------------Stalk geometry functions---------------
function [x] = xpts(N,theta,dmaj,asymmetry)
    x = zeros(1,N);
    for i = 1:N
        x(i) = (dmaj/2)*(cos(theta(i))) + asymmetry(i);
    end
end

function [y] = ypts(N,theta,dmin,asymmetry)
    y = zeros(1,N);
    for i = 1:N
        y(i) = (dmin/2)*(sin(theta(i))) + asymmetry(i);
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