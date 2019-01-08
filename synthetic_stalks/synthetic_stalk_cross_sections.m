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
%   - probably should had breadth of x and y asymmetry features.

theta = linspace(0,2*pi,100);
nloc = pi;

for i = 1:10
    
    % Basic Geometry
    min_diam = unifrnd(15,20);
    maj_diam = unifrnd(20,25);
    
    % Notch Geometry
    ndepth = unifrnd(0.1,0.25);
    nwidth = 1;
   
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
    
    % y-Asymmetryg
    aAmp = unifrnd(-0.05,0.05);
    aSym = unifrnd(-pi,pi);
    asymmetry = aAmp*sin(theta - aSym);
    
    y = min_diam*(sin(theta) + noisey + asymmetry);
    rindx = unifrnd(0.92,0.94);
    rindy = unifrnd(0.92,0.94);
    
    X(i,:) = x;
    Y(i,:) = y;
    
    plot(x,y + (i-1)*40,'k')
    hold on
    plot(x*rindx,y*rindy + (i-1)*40,'k')
    %axis([-25 25 -25 25])
    axis equal
    pause
    hold off
    
end

[coeffx,scorex,latentx,tsquaredx,explainedx,mux] = pca(X);
[coeffy,scorey,latenty,tsquaredy,explainedy,muy] = pca(Y);

subplot(1,2,1), bar(explainedx), ylim([0 100])
subplot(1,2,2), bar(explainedy), ylim([0 100])

figure
subplot(1,2,1)
plot(theta, cos(theta),'k')
hold on
plot(theta, coeffx(:,1))
plot(theta, coeffx(:,2))
plot(theta, coeffx(:,3))

subplot(1,2,2)
plot(theta, sin(theta),'k')
hold on
plot(theta, coeffy(:,1))
plot(theta, coeffy(:,2))
plot(theta, coeffy(:,3))


figure;
for i = 0:10 
    plot(coeffx(:,1) + i/50*coeffx(:,3), coeffy(:,1))
    title('PC3 - notch')
    %axis([-0.1 0.2 -0.1, 0.05])
    pause(0.25)
end



