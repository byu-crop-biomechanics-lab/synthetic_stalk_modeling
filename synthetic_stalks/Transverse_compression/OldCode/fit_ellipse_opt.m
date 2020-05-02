function [xopt, fopt, exitflag, output] = fit_ellipse_opt(Tdata,Rdata)
    % Using optimization routine (fmincon) to get the best fit of the
    % polar parametric ellipse equations to real data (with notch region
    % removed from data vectors). Assumes the ellipse has already been
    % rotated so the major axis lies upon the x axis
    
    
    % ------------Starting point and bounds------------
    %var= dmaj  dmin        % Design variables
    x0 = [100,   50];       % Starting point
    ub = [400,   300];      % Upper bound
    lb = [0.001, 0.001];    % Lower bound

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        % Design variables
        dmaj =  x(1);
        dmin =  x(2);
        
        % Other analysis variables
%         Npts = length(Rdata);
        
        % Analysis functions
        r = rpts_ellipse(Tdata,dmaj,dmin);
        
        % Objective function
        f = getfitpolar(Rdata,r); % Minimize the fit metric
        
        % Inequality constraints
        c = dmin - dmaj;
        
        % Equality constraints
        ceq = [];
        
    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'display', 'off');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
%     xopt
%     fopt
    [~,c,~] = objcon(xopt);
%     c
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end

    % ------------Stalk geometry functions---------------
%     function [x] = xpts(N,theta,dmaj)
%         x = zeros(1,N);
%         for i = 1:N
%             x(i) = (dmaj/2)*(cos(theta(i)));
%         end
%     end
% 
%     function [y] = ypts(N,theta,dmin)
%         y = zeros(1,N);
%         for i = 1:N
%             y(i) = (dmin/2)*(sin(theta(i)));
%         end
%     end
% 
%     function [xrotate,yrotate] = rotate(x,y,psi,N)
%         xrotate = zeros(1,N);
%         yrotate = zeros(1,N);
% 
%         R = [cos(psi) -sin(psi); sin(psi) cos(psi)];
% 
%         for i = 1:N
%             temp = [x(i);y(i)];
%             temp = R*temp;
%             xrotate(i) = temp(1);
%             yrotate(i) = temp(2);
%         end
%     end
    
    
%     function [fitmetric] = getfit(xreal,xsynth,yreal,ysynth,zreal,zsynth)
%         % Get the overall fit of a synthetic shape against a real shape, in a
%         % single value similar to a standard deviation
% 
%         N = length(xreal);
% 
%         if nargin == 2
%             yreal = zeros(N);
%             ysynth = zeros(N);
%             zreal = zeros(N);
%             zsynth = zeros(N);
% 
%         elseif nargin == 4
%             zreal = zeros(N);
%             zsynth = zeros(N);
% 
%         elseif mod(nargin,2) == 1
%             error('Odd number of arguments');
%         end
% 
%         sum = 0;
% 
%         for i = 1:N
%             xsq = (xsynth(i) - xreal(i))^2;
%             ysq = (ysynth(i) - yreal(i))^2;
%             zsq = (zsynth(i) - zreal(i))^2;
%             sum = sum + xsq + ysq + zsq;
%         end
% 
%         fitmetric = sqrt(sum/N);
% 
%     end
    
    
    function [r] = rpts_ellipse(theta,dmaj,dmin)
        N = length(theta);
        r = zeros(1,N);
        for i = 1:N
            r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
                + ((dmaj/2)*sin(theta(i)))^2);
        end
    end
    
    
    function [fitmetric] = getfitpolar(Rreal,rsynth)
        % Get the overall fit of a synthetic shape against a real shape, in a
        % single value similar to a standard deviation

        N = length(Rreal);

        sum = 0;

        for i = 1:N
            rsq = (rsynth(i) - Rreal(i))^2;
            sum = sum + rsq;
        end

        fitmetric = sqrt(sum/N);

    end



end
