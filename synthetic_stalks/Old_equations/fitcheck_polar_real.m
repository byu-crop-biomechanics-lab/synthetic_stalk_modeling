% Visually check the fit produced by stalk_cross_fit.m function
load AC_SP.mat ext_T ext_R
% numsection = 1;
Tdata = ext_T;
Rdata = ext_R;
N = length(Rdata);

[xopt, fopt, exitflag, output] = stalk_cross_fit_polar(Tdata,Rdata);

% Get variable values from the output of stalk_cross_fit.m:
dmaj =      xopt(1);
dmin =      xopt(2);
ndepth =    xopt(3);
nwidth =    xopt(4);
nloc =      xopt(5);
aAmp =      xopt(6);
aSym =      xopt(7);

% Analysis functions
asymmetry = aAmp*sin(Tdata - aSym);
notch = notch_fn(N,ndepth,nwidth,nloc,Tdata);
rsynth = rpts(N,Tdata,dmaj,dmin,asymmetry,notch);

% Plot the original data and the fit data on top of each other
polarplot(Tdata,Rdata);
hold on
polarplot(Tdata,rsynth)
hold off
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