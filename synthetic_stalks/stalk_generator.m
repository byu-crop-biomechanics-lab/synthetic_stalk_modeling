% Ryan Larson - 3/14/2019
%
% Reformulation of GenerateStalkSegmentBetsy_v3.m using function methods,
% corrected equations for stalk cross section shape, and outputs data in an
% array that will be easy to work with for PCA of the longitudinal data.

clear
clc

%% Create data structure for holding lots of stalk information
N = 100;    % Number of points around the stalk cross section
n = 5;    % Number of stalks to generate
nslices = 120;  % Number of z slices per stalk
theta = linspace(0,2*pi,N);

aAmplim = 0.05;

% Use radius-z scheme for handling data so it doesn't get huge
STALKS = zeros(n,nslices,N);

for numstalk = 1:n
    % Generate realish looking stalk with synthetic data
    %
    % Features include 3 segments, upper segments taper, bowed segments,
    % full interior, and rotated middle section.
    %
    % JH 8/8/2018
    %
    % =================== FROM DR COOK's "synthetic_stalk_sections.m" =====
    % parameterized corn stalk cross-section
    % Enumerating parameters:
    %   - major diameter
    %   - minor diameter
    %   - notch depth
    %   - notch width
    %   - amplitude of x-asymmetry
    %   - location of x-asymmetry
    %   - amplitude of y-asymmetry
    %   - location of y-asymmetry
    %   - probably should add breadth of x and y asymmetry features.

    clear
    clc

    theta = linspace(0,2*pi,100);
    nloc = pi;
    m = 1; % multiplier

    % Data = [1;
    %     24;
    %     22.53;
    %     23.3;
    %     21;
    %     0;
    %     0;
    %     41.08];

    for i = 1:3

        %     % Basic Geometry
        %     min_diam = unifrnd(15,20);
        %     maj_diam = unifrnd(20,25);
        %    % Notch Geometry
        %     ndepth = unifrnd(0.1,0.25);
        %     nwidth = 1;

    %     if i == 1     % Segments 1-3 (Bottom)
    %         min_diam = 22.53;
    %         maj_diam = 24;
    %         ndepth = 0;
    %         nwidth = 1;
    %     elseif i == 2
    %         min_diam = 24.84;
    %         maj_diam = 26;
    %         ndepth = .026;
    %         nwidth = 1;
    %     else
    %         min_diam = 27.11;
    %         maj_diam = 29.51;
    %         ndepth = .025;
    %         nwidth = 1;
    %     end
    %     

    % DEFINE DESIRED GEOMETRY FOR EACH INTERNODE
            if i == 1
                min_diam = 22.6;
                maj_diam = 26.84;
                ndepth = 0.023;
                nwidth = 1;
            elseif i == 2
                min_diam = 21.3;
                maj_diam = 23.75;
                ndepth = .168;
                nwidth = 1;
            else
                min_diam = 19.64;
                maj_diam = 21.15;
                ndepth = .303;
                nwidth = 1;
            end



        min_diam = min_diam*m/2;
        maj_diam = maj_diam*m/2;
        ndepth = ndepth*m;


        % x-Asymmetry
        aAmp = unifrnd(-0.05,0.05);
        aSym = unifrnd(-pi,pi);

        % Random Noise
        noisex = unifrnd(-0.001,0.001,1,100);
        noisey = unifrnd(-0.001,0.001,1,100);

        notch = ndepth./cosh(4*(theta/nwidth-nloc)).^2;
        asymmetry = aAmp*sin(theta - aSym);

        x = maj_diam*(cos(theta) + noisex + notch + asymmetry);

        %if mod(i,2) == 0, x = x*-1; end

        % y-Asymmetry
        aAmp = unifrnd(-0.05,0.05);
        aSym = unifrnd(-pi,pi);
        asymmetry = aAmp*sin(theta - aSym);

        y = min_diam*(sin(theta) + noisey + asymmetry);
        rindx = unifrnd(0.92,0.94);
        rindy = unifrnd(0.92,0.94);

        X(i,:) = x;
        Y(i,:) = y;
    end

    %========================================================================
    alpha = -(2*pi)*.5;         % must be a negative number for the shift function
    % to work properly
    % rotation matrix reassignment
    new = ([X(2,:)' Y(2,:)']) *([cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]);
    Y(2,:) = new(:,2)';
    X(2,:) = new(:,1)';


    % Shift function
    % Reassigns index values to compensate for the rotation
    shift = round((-100*alpha)/(2*pi));

    X2prime = zeros(1,(length(X) + shift));
    X2prime((shift+1):end) = X(2,:);
    X2prime(1:shift+1) = X2prime((length(X):end));
    X2prime = X2prime(1:length(X));
    X(2,:) = X2prime;

    Y2prime = zeros(1,(length(Y) + shift));
    Y2prime((shift+1):end) = Y(2,:);
    Y2prime(1:shift+1) = Y2prime((length(Y):end));
    Y2prime = Y2prime(1:length(Y));
    Y(2,:) = Y2prime;

    r = 15;
    % Plot cross sections of stalks
    hold on
    set(gcf, 'Position', [300, 300, 1000, 350])

    subplot(1,3,1)
    plot(X(1,:),Y(1,:))
    xlim([-r r])
    ylim([-r r])

    subplot(1,3,2)
    plot(X(2,:),Y(2,:))
    xlim([-r r])
    ylim([-r r])

    subplot(1,3,3)
    plot(X(3,:),Y(3,:))
    xlim([-r r])
    ylim([-r r])

    Xs(1,:) = X(1,:);
    Ys(1,:) = Y(1,:);
    Xs(1,100) = Xs(1,1);    % Solve the small gap problem
    Ys(1,100) = Ys(1,1);

    Xs(2,:) = X(2,:);
    Ys(2,:) = Y(2,:);
    Xs(2,100) = Xs(2,1);
    Ys(2,100) = Ys(2,1);

    Xs(3,:) = X(3,:);
    Ys(3,:) = Y(3,:);
    Xs(3,100) = Xs(3,1);
    Ys(3,100) = Ys(3,1);

    %use for loop to create multiple versions with different internode lengths
    % for q=1:2

    % Populate matrices (this method relies on the fact that the X and Y coords
    % are centered at the origin)
    L = 83*m; % internode length
    a = L/40; % multiplier

    Z = zeros(120,100);

    for i = 1:120       % Generate Z vector
        Z(i,:) = i*a; %3*q instead of 4 if you are varying the lengths
    end

    % First segment
    for i = 1:40
        X(i,:) = (.00000125*(i-20)^4 +.95)*Xs(1,:);
        Y(i,:) = (.00000125*(i-20)^4 +.95)*Ys(1,:);
    end

    % Second Segment

    for i = 41:80
        X(i,:) = (.00000125*(i-60)^4 +.9)*Xs(2,:);
        Y(i,:) = (.00000125*(i-60)^4 +.9)*Ys(2,:);
    end


    % Third Segment
    
    for i = 81:120
        X(i,:) = (.00000125*(i-100)^4 +.85)*Xs(3,:);
        Y(i,:) = (.00000125*(i-100)^4 +.85)*Ys(3,:);
    end


    surf2stl_V1('stalk.stl',X,Y,Z)
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
