clear all;
close all;

% syms theta min_diam maj_diam nloc ndepth nwidth x y
N = 1000;
theta = linspace(0,2*pi,N);
min_diam = 22;
maj_diam = 27;
nloc = pi;
ndepth = linspace(0.001,1,100);
nwidth = linspace(0.001,1,100);


notch = zeros(length(ndepth),length(nwidth),length(theta));
x = zeros(length(ndepth),length(nwidth),length(theta));
y = zeros(length(ndepth),length(nwidth),length(theta));
dxdy = zeros(length(ndepth),length(nwidth),length(theta));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        for k = 1:length(theta)
            notch(i,j,k) = ndepth(i)/cosh(4*(theta(k)/nwidth(j)-nloc))^2;
            x(i,j,k) = maj_diam*(cos(theta(k)) + notch(i,j,k));
            y(i,j,k) = min_diam*(sin(theta(k)));
            dxdy(i,j,k) = -(maj_diam*(sin(theta(k)) - (8.*ndepth(i)*sinh(4*nloc - ...
                (4*theta(k))/nwidth(j)))/(nwidth(j)*cosh(4*nloc - ...
                (4*theta(k))/nwidth(j))^3)))/(min_diam*cos(theta(k)));
        end
    end
end

dxdy(:,:,1) = 0;

% x = maj_diam.*(cos(theta) + notch);
% y = min_diam.*(sin(theta));

% dxdy = -(maj_diam.*(sin(theta) - (8.*ndepth.*sinh(4.*nloc - (4.*theta)./nwidth))./(nwidth.*cosh(4.*nloc - (4.*theta)./nwidth).^3)))./(min_diam.*cos(theta));
% dxdy(1) = 0;

% Calculate the theta values where the dxdy slope is zero (measurement
% points for nwidth and ndepth)
thetapts = zeros(length(ndepth),length(nwidth),3);
% thetaptstemp = zeros(3,1);
dxdytemp = zeros(1,length(theta));
count = 1;
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        for k = 1:(length(theta)-2)
            % If there's a sign change, catch it as a possible measurement
            % location
            if (dxdy(i,j,k)<0 & dxdy(i,j,k+1)>0) || (dxdy(i,j,k)>0 & dxdy(i,j,k+1)<0)
                % If the possible location is not an asymptote, catch it
                if abs(dxdy(i,j,k)-dxdy(i,j,k+1)) < 10
                    % Create a temporary holding vector for dxdy across all
                    % theta, for each specific combination of nwidth and
                    % ndepth
                    for m = 1:(N-2)
                        dxdytemp(m) = dxdy(i,j,m);
                    end
%                     dxdytemp = dxdy(i,j,:);

                    % Fill a temporary vector of theta values for
                    % measurement points
                    thetaptstemp = interp1(dxdytemp((k-2):(k+3)),theta((k-2):(k+3)),0,'pchip');
                    thetapts(i,j,count) = thetaptstemp;
                    
%                     thetapts(i,j,count) = interp1(dxdytemp(k-4:k+5),theta(k-4:k+5),0,'pchip');
                    count = count + 1;
                    if count > 3
                        count = 1;
                    end
%                     thetapts(i,j,count) = interp1(dxdy(i,j,k-4:k+5),theta(i,j,k-4:k+5),0,'pchip');
                    
                    
%                     theta_zero = interp1(dxdy((i,j,k-4):(i,j,k+5)),theta((i,j,k-4):(i,j,k+5)),0,'pchip');
%                     thetapts = [Zeros;theta_zero];
                end
            end
        end
    end
end

% Establish x and y values for measurement points
xmeasurepts = zeros(size(thetapts));
ymeasurepts = zeros(size(thetapts));
notchmeasurepts = zeros(size(thetapts));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        for k = 1:3
            notchmeasurepts(i,j,k) = ndepth(i)/cosh(4*(thetapts(i,j,k)/nwidth(j)-nloc)).^2;
            xmeasurepts(i,j,k) = maj_diam*(cos(thetapts(i,j,k)) + notchmeasurepts(i,j,k));
            ymeasurepts(i,j,k) = min_diam*(sin(thetapts(i,j,k)));
        end
    end
end
        
% for i = 1:length(zeros)
%     notchmeasurepts(i) = ndepth./cosh(4*(zeros(i)/nwidth-nloc)).^2;
%     xmeasurepts(i) = maj_diam*(cos(zeros(i)) + notchmeasurepts(i));
%     ymeasurepts(i) = min_diam*(sin(zeros(i)));
% end

% Measure the true notch depth and width
ndepth_true = zeros(length(ndepth),length(nwidth));
nwidth_true = zeros(length(ndepth),length(nwidth));
for i = 1:length(ndepth)
    for j = 1:length(nwidth)
        ndepth_true(i,j) = (abs(xmeasurepts(i,j,1)-xmeasurepts(i,j,2)) + abs(xmeasurepts(i,j,2)-xmeasurepts(i,j,3)))/2;
        nwidth_true(i,j) = abs(ymeasurepts(i,j,1)-ymeasurepts(i,j,3));
    end
end




% ndepth_true = (abs(xmeasurepts(1)-xmeasurepts(2)) + abs(xmeasurepts(2)-xmeasurepts(3)))/2
% nwidth_true = abs(ymeasurepts(1)-ymeasurepts(3))

% figure(1)
% plot(theta,dxdy)
% 
% figure(2)
% plot(x,y)


%% OLD CODE
% dmaj = 2;
% dmin = 1;
% 
% nwidth = 0.5;
% ndepth = 0.1;
% if nwidth >= dmaj
%     error('Notch width is too large. Increase dmaj or decrease nwidth.')
% end
% 
% theta = linspace(-pi,pi,500);   % The range of theta must stay as -pi to pi
% yshift = zeros(1,length(theta));
% notch = zeros(1,length(theta));
% y = (dmin/2)*cos(theta);
% x = (dmaj/2)*sin(theta);
% % flatlimit = pi/4;
% flatlimit = asin(nwidth/dmaj);
% for i = 1:length(notch)
%     if theta(i) < -flatlimit | theta(i) > flatlimit
%         notch(i) = 0;
%     else
%         notch(i) = ((ndepth/2)*(cos((2*pi/nwidth)*theta(i)))+(ndepth/4));
%     end
% end
% 
% % flatheight = 0;
% 
% for i = 1:length(yshift)
%     % Don't shift y if outside the chosen bounds
%     if theta(i) < -pi | theta(i) > pi
%         yshift(i) = 0;
%     elseif theta(i) >= -flatlimit & theta(i) <= flatlimit
% %         yshift(i) = (dmin/2)*cos(-flatlimit);   % Keep the same value as the value at flatlimits
%         yshift(i) = y(i) - notch(i);
%     else
%         yshift(i) = y(i);
%     end          
% end
% 
% % yshift(yshift == 0) = NaN;  % THIS STILL CAUSES A HOLE IN THE BOTTOM OF THE ELLIPSE
% % for i = 1:length(theta)
% %     if theta(i) == 0
% %         yshift(i) = 0;
% %     end
% % end
% % x(1) = 0;
% % x(end) = 0;
% 
% % figure(1)
% % plot(theta,y)
% % hold on
% % plot(theta,yshift)
% % legend('y','yshift');
% % xlabel('\theta, radians')
% 
% hold off
% figure(2)
% plot(x,yshift)
% 
% figure(3)
% plot(theta,notch)