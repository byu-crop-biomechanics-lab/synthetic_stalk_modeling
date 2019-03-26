% Ryan Larson - 3/14/2019
%
% Reformulation of GenerateStalkSegmentBetsy_v3.m using function methods,
% corrected equations for stalk cross section shape, and outputs data in an
% array that will be easy to work with for PCA of the longitudinal data.

clear
clc

%% Create data structure for holding lots of stalk information
N = 100;    % Number of points around the stalk cross section
n = 300;    % Number of stalks to generate
nslices = 120;  % Number of z slices per stalk
theta = linspace(0,2*pi,N);

aAmplim = 0.05;

% Use radius-z scheme for handling data so it doesn't get huge
STALKS = zeros(n,nslices,N);

for numstalk = 1:n
    


    %% Generate stalks
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



    nloc = pi;
    m = 1; % multiplier

    for i = 1:3

        %     % Basic Geometry
        %     min_diam = unifrnd(15,20);
        %     maj_diam = unifrnd(20,25);
        %    % Notch Geometry
        %     ndepth = unifrnd(0.1,0.25);
        %     nwidth = 1;


    % DEFINE DESIRED GEOMETRY FOR EACH INTERNODE
            if i == 1
                min_diam = unifrnd(19,20);
                maj_diam = unifrnd(23,25);
                ndepth = unifrnd(0.1,0.25);
                nwidth = unifrnd(1,4);
            elseif i == 2
                min_diam = unifrnd(17,19);
                maj_diam = unifrnd(22,23);
                ndepth = unifrnd(0.1,0.25);
                nwidth = unifrnd(1,4);
            else
                min_diam = unifrnd(15,17);
                maj_diam = unifrnd(20,22);
                ndepth = unifrnd(0.1,0.25);
                nwidth = unifrnd(1,4);
            end



        min_diam = min_diam*m/2;
        maj_diam = maj_diam*m/2;
        ndepth = ndepth*m;


        % Random Noise
        noisex = unifrnd(-0.005,0.005,1,N);
        noisey = unifrnd(-0.005,0.005,1,N);

        notch = notch_fn(N,ndepth,nwidth,nloc,theta);
        asymmetry = Asymmetry(aAmplim,theta,N);

        x = xpts(N,theta,notch,maj_diam,noisex,asymmetry);

        %if mod(i,2) == 0, x = x*-1; end

        % y-Asymmetry
        aAmp = unifrnd(-0.05,0.05);
        aSym = unifrnd(-pi,pi);
        asymmetry = aAmp*sin(theta - aSym);

        y = ypts(N,theta,min_diam,noisey,asymmetry);
        
        X = zeros(3,N);
        Y = zeros(3,N);
        
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
    shift = round((-N*alpha)/(2*pi));

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
        X(i,:) = (.00000125*(i-20)^4 + 1)*Xs(1,:);
        Y(i,:) = (.00000125*(i-20)^4 + 1)*Ys(1,:);
    end

    % Second Segment
    for i = 41:80
        X(i,:) = (.00000125*(i-60)^4 + 1)*Xs(2,:);
        Y(i,:) = (.00000125*(i-60)^4 + 1)*Ys(2,:);
    end

    % Third Segment
    for i = 81:120
        X(i,:) = (.00000125*(i-100)^4 + 1)*Xs(3,:);
        Y(i,:) = (.00000125*(i-100)^4 + 1)*Ys(3,:);
    end

    % % This block is necessary because of the nature of the STL writer.
    % % In the stitching process, it will make triangle connections from one
    % % index to the other.  The first 60 indices are responsible for the
    % % exterior where the generator works upwards.  The next 60 indices must
    % % work from 'top to bottom'
    % flipper = flip(Z);
    % Z(61:120,:) = flipper;
    % 
    % X(61:120,:) = zeros;
    % flipper = flip(X(1:60,:));
    % X(61:120,:) = .8*flipper;
    % 
    % Y(61:120,:) = zeros;
    % flipper = flip(Y(1:60,:));
    % Y(61:120,:) = .8*flipper;
    % 
    % Z(121,:) = 4;
    % X(121,:) = X(1,:);
    % Y(121,:) = Y(1,:);

    % also note the necessity of duplicate points to finish the stitching

    %% Convert Cartesian data to polar
    [~,R] = cart2pol(X,Y);

    % Insert new data into STALKS
    for i = 1:N
        for j = 1:nslices
            STALKS(numstalk,j,i) = R(j,i);
        end
    end

    
    %% Write to STL if desired
    % surf2stl_V1('stalk.stl',X,Y,Z)

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
