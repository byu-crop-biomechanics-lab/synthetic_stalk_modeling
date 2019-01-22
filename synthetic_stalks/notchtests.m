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
            notch(i,j,k) = ndepth(i)./cosh((10/nwidth(j))*(theta(k)-nloc)).^2;
            x(i,j,k) = maj_diam*(cos(theta(k)) + notch(i,j,k));
            y(i,j,k) = min_diam*(sin(theta(k)));
            dxdy(i,j,k) = -(maj_diam*(sin(theta(k)) - (8.*ndepth(i)*sinh(4*nloc - ...
                (4*theta(k))/nwidth(j)))/(nwidth(j)*cosh(4*nloc - ...
                (4*theta(k))/nwidth(j))^3)))/(min_diam*cos(theta(k)));
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


thetapts = nan(length(ndepth),length(nwidth),3);
dxdytemp = nan(1,length(theta));
count = 1;
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        for k = 1:length(theta)
            % NEW SOLUTION: SEARCH IN EACH OF THE EXPECTED ANGLE RANGES FOR
            % A MEASUREMENT FEATURE. IF THERE ISN'T A MEASUREMENT FEATURE,
            % THEN LET ALL ENTRIES FOR THAT NDEPTH/NWIDTH COMBO REMAIN NAN
            
            if theta(k)>pi/2 & theta(k)<(3*pi)/2
                % If there's a sign change, catch it
                if (dxdy(i,j,k)<0 & dxdy(i,j,k+1)>0) || (dxdy(i,j,k)>0 & dxdy(i,j,k+1)<0)
                    % Create a temporary holding vector for dxdy across all
                    % theta, for each specific combination of nwidth and
                    % ndepth
                    for m = 1:(N-2)
                        dxdytemp(m) = dxdy(i,j,m);
                    end

                    % Fill a temporary vector of theta values for
                    % measurement points
                    thetaptstemp = interp1(dxdytemp((k-1):(k+1)),theta((k-1):(k+1)),0,'pchip');
                    thetapts(i,j,count) = thetaptstemp;

                    count = count + 1;
                    if count > 3
                        count = 1;
                    end
                end
            end
        end
    end
end

% % Troubleshooting: Check how many ndepth/nwidth combinations don't have all
% % three theta slots filled
% problemvecs = 0;
% for i = 1:length(ndepth)
%     for j = 1:length(nwidth)
%         if nnz(~thetapts(i,j,:)) > 0
%             problemvecs = problemvecs+1;
%         end
%     end
% end
% problemratio = problemvecs/(N^2);
            
% Establish x and y values for measurement points
xmeasurepts = nan(size(thetapts));
ymeasurepts = nan(size(thetapts));
notchmeasurepts = nan(size(thetapts));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        for k = 1:3
            notchmeasurepts(i,j,k) = ndepth(i)/cosh(4*(thetapts(i,j,k)/nwidth(j)-nloc))^2;
            xmeasurepts(i,j,k) = maj_diam*(cos(thetapts(i,j,k)) + notchmeasurepts(i,j,k));
            ymeasurepts(i,j,k) = min_diam*(sin(thetapts(i,j,k)));
        end
    end
end


% % TRY PLOTTING THE CROSS SECTION SHAPES TO SEE HOW THEY CHANGE AS A
% % FUNCTION OF NDEPTH OR NWIDTH
% %// Plot starts here
figure

% ADDITION: MAKE THE BELOW PLOT INTO A SUBPLOT SO YOU CAN COMPARE THE
% RESULTS OF DX/DY AGAINST THETA AT THE SAME TIME (CATCH ASYMPTOTES AND
% ZEROS AND MAKE SURE THEY CORRESPOND)


%// Plot point by point
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        xtemp = zeros(1,N);
        ytemp = zeros(1,N);
        for k = 1:length(theta)
            xtemp(k) = x(i,j,k);
            ytemp(k) = y(i,j,k);
        end
        plot(xtemp,ytemp)
        hold on
        scatter(xmeasurepts(i,j,:),ymeasurepts(i,j,:))
        % plot points that show the measurement points
        hold off
        pause(0.05);
    end
end

% close(figure(1))
% %// Plot starts here
% figure(2)
% 
% %// Plot point by point
% for i = 1:5
%     for j = 1:length(nwidth)
%         notchtemp = zeros(1,N);
%         for k = 1:length(theta)
%             notchtemp(k) = notch(i,j,k);
%         end
%         plot(theta,notchtemp)
%         pause(0.05);
%     end
% end
% close(figure(2))





% Measure the true notch depth and width
ndepth_true = NaN(length(ndepth),length(nwidth));
nwidth_true = NaN(length(ndepth),length(nwidth));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        ndepth_true(i,j) = ((xmeasurepts(i,j,1)-xmeasurepts(i,j,2)) + (xmeasurepts(i,j,3)-xmeasurepts(i,j,2)))/2;
        nwidth_true(i,j) = ymeasurepts(i,j,1)-ymeasurepts(i,j,3);
    end
end

% % Visualize with 3D plots
% figure(1)
% surf(ndepth_true)
% xlabel('ndepth')
% ylabel('nwidth')
% zlabel('True notch depth')
% figure(2)
% surf(nwidth_true)
% xlabel('ndepth')
% ylabel('nwidth')
% zlabel('True notch width')

% ndepth_true = (abs(xmeasurepts(1)-xmeasurepts(2)) + abs(xmeasurepts(2)-xmeasurepts(3)))/2
% nwidth_true = abs(ymeasurepts(1)-ymeasurepts(3))

% figure(1)
% plot(theta,dxdy)
% 
% figure(2)
% plot(x,y)