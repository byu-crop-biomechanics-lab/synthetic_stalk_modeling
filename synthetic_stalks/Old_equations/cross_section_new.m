% cross_section_new.m: Script to test out new way of defining the cross
% sections so that the notch can be located and sized independently of the
% major and minor diameters
clear;
close;

% Define main parameters
N = 100;
theta = linspace(0,2*pi,N);

dmaj = 40;
dmin = 15;

ndepth = 8;
nwidth = 5;

for j = 1:3
    for i = linspace(pi/2,3*pi/2,100)
        nloc = i;
        phi = nloc - pi;

        xnotch = notch_fn(N,ndepth,nwidth,nloc,theta);
        ynotch = zeros(1,N);

        [xnotch_rot,ynotch_rot] = rotate(xnotch,ynotch,phi,N);

        xmain = dmaj*cos(theta);
        ymain = dmin*sin(theta);

        x = xmain + xnotch_rot;
        y = ymain + ynotch_rot;

        plot(x,y);
        xlim([-(dmaj + 10),(dmaj + 10)]);
        ylim([-(dmin + 10),(dmin + 10)]);
    %     hold on
        pause(0.05);

    end
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