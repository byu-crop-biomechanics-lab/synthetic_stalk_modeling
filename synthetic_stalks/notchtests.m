clear all;
close all;

N = 1000;
theta = linspace(0,2*pi,N);
min_diam = 22;
maj_diam = 27;
min_maj_ratio = min_diam/maj_diam;
nloc = pi;
ndepth = linspace(0,1,100);
nwidth = linspace(0,5,100); % limits are determined by the 


notch = zeros(length(ndepth),length(nwidth),length(theta));
x = zeros(length(ndepth),length(nwidth),length(theta));
y = zeros(length(ndepth),length(nwidth),length(theta));
dxdy = zeros(length(ndepth),length(nwidth),length(theta));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        for k = 1:length(theta)
            notch(i,j,k) = ndepth(i)/cosh((10/nwidth(j))*(theta(k)-nloc))^2;
            x(i,j,k) = maj_diam*(cos(theta(k)) + notch(i,j,k));
            y(i,j,k) = min_diam*(sin(theta(k)));
            dxdy(i,j,k) =  maj_diam*(-(20*ndepth(i)*tanh((10*(theta(k)-nloc))/nwidth(j))*(sech((10*(theta(k)-nloc))/nwidth(j)))^2)/nwidth(j) - sin(theta(k)))/(min_diam*cos(theta(k)));
        end
    end
end

% dxdy(:,:,1) = 0;    % FIXED: THIS LINE IS NOT NEEDED IF THE NEXT FIX BELOW IS IMPLEMENTED

% Calculate the theta values where the dxdy slope is zero (measurement
% points for nwidth and ndepth)

% NOTE 1/18/2019: THE FOLLOWING BLOCK OF CODE DOES NOT CATCH ALL THREE
% MEASUREMENT POINTS NEARLY HALF THE TIME (PROBLEMRATIO RETURNS 46.7%).
% THIS INDICATES A PROBLEM WITH THE METHOD. IT DOESN'T TAKE INTO ACCOUNT
% CASES WHERE THE NOTCH DEPTH IS NOT DEEP ENOUGH TO CAUSE 


nopts = 0;
yespts = 0;
thetapts = nan(length(ndepth),length(nwidth),3);
dxdytemp = nan(1,length(theta));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        
        count = 1;
        % NEW SOLUTION: SEARCH IN EACH OF THE EXPECTED ANGLE RANGES FOR
        % A MEASUREMENT FEATURE. IF THERE ISN'T A MEASUREMENT FEATURE,
        % THEN LET ALL ENTRIES FOR THAT NDEPTH/NWIDTH COMBO REMAIN NAN
        
        % Create a temporary holding vector for dxdy across all
        % theta, for each specific combination of nwidth and
        % ndepth
        for m = 1:length(theta)
            dxdytemp(m) = dxdy(i,j,m);
        end

        %% Catch the first possible measurement point:
        bound1 = pi/2;
        bound2 = pi;
        [maxind,minind] = boundinds(bound1,bound2,theta);

        % Check the bounds of the range in question to see if dx/dy
        % is NaN there
        if isnan(dxdytemp(maxind)) & isnan(dxdytemp(minind))
            thetapts(i,j,count) = NaN;
        else
            thetaptstemp = interp1(dxdytemp(minind:maxind),theta(minind:maxind),0,'pchip');

%             % Verify that the interpolated value is good for the range
%             % of interest:
%             if (thetaptstemp<theta(maxind)) & (thetaptstemp>theta(minind))
%                 thetapts(i,j,count) = thetaptstemp;
%                 count = count + 1;
%             else
%                 thetapts(i,j,count) = NaN;
%             end
            thetapts(i,j,count) = thetaptstemp;
%             disp(thetapts(i,j,count));
            count = count + 1;
        end

        %% Catch the second possible measurement point
        bound1 = pi - 0.2;
        bound2 = pi + 0.2;
        [maxind,minind] = boundinds(bound1,bound2,theta);
        
        % Check the bounds of the range in question to see if dx/dy
        % is NaN there
        if isnan(dxdytemp(maxind)) & isnan(dxdytemp(minind))
            thetapts(i,j,count) = NaN;
        else
            thetaptstemp = interp1(dxdytemp(minind:maxind),theta(minind:maxind),0,'pchip');


% %             Verify that the interpolated value is good for the range
% %             of interest:
%             if (thetaptstemp<theta(maxind)) & (thetaptstemp>theta(minind))
%                 thetapts(i,j,count) = thetaptstemp;
%                 count = count + 1;
%             else
%                 thetapts(i,j,count) = NaN;
%                 continue
%             end
            thetapts(i,j,count) = thetaptstemp;
            count = count + 1;
        end
        
        %% Catch the third possible measurement point
        bound1 = pi;
        bound2 = 3*pi/2;
        [maxind,minind] = boundinds(bound1,bound2,theta);

        % Check the bounds of the range in question to see if dx/dy
        % is NaN there
        if isnan(dxdytemp(maxind)) & isnan(dxdytemp(minind))
            thetapts(i,j,count) = NaN;
        else
            thetaptstemp = interp1(dxdytemp(minind:maxind),theta(minind:maxind),0,'pchip');


%             % Verify that the interpolated value is good for the range
%             % of interest:
%             if (thetaptstemp<theta(maxind)) & (thetaptstemp>theta(minind))
%                 thetapts(i,j,count) = thetaptstemp;
%                 count = count + 1;
%             else
%                 thetapts(i,j,count) = NaN;
%             end
            thetapts(i,j,count) = thetaptstemp;
            count = count + 1;
        end
%         disp(thetapts(i,j,:))
    end
end
            
% Establish x and y values for measurement points
xmeasurepts = nan(size(thetapts));
ymeasurepts = nan(size(thetapts));
notchmeasurepts = nan(size(thetapts));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        for k = 1:3
            notchmeasurepts(i,j,k) = ndepth(i)/cosh((10/nwidth(j))*(thetapts(i,j,k)-nloc))^2;
            xmeasurepts(i,j,k) = maj_diam*(cos(thetapts(i,j,k)) + notchmeasurepts(i,j,k));
            ymeasurepts(i,j,k) = min_diam*(sin(thetapts(i,j,k)));
        end
    end
end

% Measure the true notch depth and width
ndepth_true = NaN(length(ndepth),length(nwidth));
nwidth_true = NaN(length(ndepth),length(nwidth));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        ndepth_true(i,j) = abs((abs(xmeasurepts(i,j,1)-xmeasurepts(i,j,2)) + abs(xmeasurepts(i,j,3)-xmeasurepts(i,j,2)))/2);
        nwidth_true(i,j) = ymeasurepts(i,j,1)-ymeasurepts(i,j,3);
        disp(ndepth_true(i,j));
        disp(nwidth_true(i,j));
    end
end

% % % TRY PLOTTING THE CROSS SECTION SHAPES TO SEE HOW THEY CHANGE AS A
% % % FUNCTION OF NDEPTH OR NWIDTH
% % %// Plot starts here
% figure
% thetaptsscatter = zeros(3,1);
% dxdyscatter = zeros(3,1);
% %// Plot point by point
% for i = 1:length(ndepth)
%     for j = 1:length(nwidth)
%         xtemp = zeros(1,N);
%         ytemp = zeros(1,N);
%         dxdytemp = zeros(1,N);
%         for k = 1:length(theta)
%             xtemp(k) = x(i,j,k);
%             ytemp(k) = y(i,j,k);
%             dxdytemp(k) = dxdy(i,j,k);
%         end
%         for k = 1:3
%             thetaptsscatter(k) = thetapts(i,j,k);
%             dxdyscatter(k) = maj_diam*(-(20*ndepth(i)*tanh((10*(thetaptsscatter(k)-nloc))...
%                 /nwidth(j))*(sech((10*(thetaptsscatter(k)-nloc))/nwidth(j)))^2)/nwidth(j) - ...
%                 sin(thetaptsscatter(k)))/(min_diam*cos(thetaptsscatter(k)));
%         end
%         
%         
%         subplot(1,2,1);
%         plot(xtemp,ytemp)
%         hold on
%         scatter(xmeasurepts(i,j,:),ymeasurepts(i,j,:));  % plot measurement points over the top of the cross section shape
%         hold off
%         title('Cross sectional shape and measurement points');
%         subplot(1,2,2);
%         plot(theta(270:730),dxdytemp(270:730))
%         hold on
%         scatter(thetaptsscatter,dxdyscatter);
%         hold off
%         xlabel('\theta, rad');
%         ylabel('dx/dy')
%         pause(0.05);
%     end
% end


% Visualize with 3D plots
figure(1)
surf(ndepth_true)
xlabel('ndepth')
ylabel('nwidth')
zlabel('True notch depth')
figure(2)
surf(nwidth_true)
xlabel('ndepth')
ylabel('nwidth')
zlabel('True notch width')

% ndepth_true = (abs(xmeasurepts(1)-xmeasurepts(2)) + abs(xmeasurepts(2)-xmeasurepts(3)))/2
% nwidth_true = abs(ymeasurepts(1)-ymeasurepts(3))

figure(1)
plot(theta,dxdy)

figure(2)
plot(x,y)