function [R_int] = normint(R,theta,t)
    % Calculates the new interior points based on local central difference
    % derivative
    close;
    % Convert to Cartesian
    [xy_columns] = convert_to_xy(R,theta);
    x0 = xy_columns(:,1);
    y0 = xy_columns(:,2);
    figure(1);
    plot(x0,y0,'b.');
    hold on
    plot(x0(1),y0(1),'r.');
    
    % Get vector of normal slopes
    slopes = zeros(size(theta));
    m = zeros(size(theta));
    
    for i = 1:length(slopes)
        if i == 1
            m(i) = (y0(2)-y0(end))/(x0(2)-x0(end));
        elseif i == length(slopes)
            m(i) = (y0(1)-y0(end-1))/(x0(1)-x0(end-1));            
        else
            m(i) = (y0(i+1)-y0(i-1))/(x0(i+1)-x0(i-1));
        end
        
    end
    
    for i = 1:length(slopes)
        slopes(i) = -1/m(i);
    end
    
    % Compute normal slope angles
    phi = zeros(size(slopes));
    for i = 1:length(phi)
        phi(i) = atan(slopes(i));
    end
    
    % Compute locations of interior points along the normal vectors
    x1 = zeros(size(x0));
    y1 = zeros(size(y0));
    rise = zeros(size(y0));
    run = zeros(size(x0));
    
    for i = 1:length(slopes)
        run(i) = t*cos(phi(i));
        rise(i) = t*sin(phi(i));
        xtemp1 = x0(i) - run(i);
        xtemp2 = x0(i) + run(i);
        ytemp1 = y0(i) - rise(i);
        ytemp2 = y0(i) + rise(i);
        if abs(xtemp1) <= abs(xtemp2)
            x1(i) = xtemp1;
        else
            x1(i) = xtemp2;
        end
        if abs(ytemp1) <= abs(ytemp2)
            y1(i) = ytemp1;
        else
            y1(i) = ytemp2;
        end
%         y1(i) = y0(i) - rise(i);
    end
    
    plot(x1,y1,'r.');
    for i = 1:length(x1)
        plot([x0(i) x1(i)],[y0(i) y1(i)],'g');
    end
    axis equal
    hold off
    
    % Convert back to polar coordinates
    [thetatemp, Rtemp] = cart2pol(x1,y1);
    thetatemp = wrapTo2Pi(thetatemp);   % Make all the negative pi values sit on a 0-2*pi system
    
    R_int_interp = interp1(thetatemp,Rtemp,theta,'pchip','extrap'); 
    R_int = R_int_interp;
    figure(2);
    polarplot(theta,R,'b.');
    hold on
    polarplot(thetatemp,R_int,'r.');
    hold off

end


% function [r] = rpts(N,theta,dmaj,dmin)
%     % theta is in radians
%     r = zeros(1,N);
%     for i = 1:N
%         r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
%             + ((dmaj/2)*sin(theta(i)))^2);
%     end
% end


function [xy_columns] = convert_to_xy(R,theta)
    N = length(theta);
    xy_columns = zeros(N,2);
    for i = 1:N
        xy_columns(i,1) = R(i)*cos(theta(i));
        xy_columns(i,2) = R(i)*sin(theta(i));
    end
end