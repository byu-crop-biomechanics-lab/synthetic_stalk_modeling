function [intpts] = normint(R,theta,t)
    % Calculates the new interior points based on local central difference
    % derivative
    
    % Convert to Cartesian
    [xy_columns] = convert_to_xy(R,theta);
    x = xy_columns(:,1);
    y = xy_columns(:,2);
    
    % Get vector of normal slopes
    slopes = zeros(size(theta));
    
    for i = 1:length(slopes)
        if i == 1
            slopes(i) = (y(2)-y(end))/(x(2)-x(end));
        elseif i == length(slopes)
            slopes(i) = (y(1)-y(end-1))/(x(1)-x(end-1));            
        else
            slopes(i) = (y(i+1)-y(i-1))/(x(i+1)-x(i-1));
        end
        
    end



end


function [r] = rpts(N,theta,dmaj,dmin)
    % theta is in radians
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end


function [xy_columns] = convert_to_xy(R,theta)
    N = length(theta);
    xy_columns = zeros(N,2);
    for i = 1:N
        xy_columns(i,1) = R(i)*cos(theta(i));
        xy_columns(i,2) = R(i)*sin(theta(i));
    end
end