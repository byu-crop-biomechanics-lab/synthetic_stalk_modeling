function [xopt, fopt, exitflag, output] = stalk_cross_fit_polar(numsection)
    % Using optimization routine (fmincon) to get the best fit of the
    % synthetic stalk equations to an example shape
    
%     load polar_cross_sections.mat sections
    load TRryan.mat rhoDCSR
    
%     R = sections(numsection,:);
    R = rhoDCSR(1,:,numsection);
    
    % ------------Starting point and bounds------------
    %var= dmaj  dmin  ndepth  nwidth  nloc    aAmp   aSym  % Design variables
    x0 = [20,   10,   5,      1.1,    pi,     0,     0];   % Starting point
    ub = [30,   20,   10,     5,      3*pi/2, 10,  pi];  % Upper bound
    lb = [0.001,0.001,0.01,   0.01,   pi/2,   -10, -pi]; % Lower bound

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
        ndepth =  x(3);
        nwidth = x(4);
        nloc = x(5);
        aAmp = x(6);
        aSym = x(7);
        
        % Other analysis variables
        N = length(R);
        theta = linspace(0,2*pi,N);
        
        % Analysis functions
        asymmetry = aAmp*sin(theta - aSym);
        notch = notch_fn(N,ndepth,nwidth,nloc,theta);
        r = rpts(N,theta,dmaj,dmin,asymmetry,notch);
        
        % Objective function
        f = getfitpolar(R,r); % Minimize the fit metric
        
        % Inequality constraints
        c = dmin - dmaj;
        
        % Equality constraints
        ceq = [];
        
    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'display', 'iter-detailed');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    xopt
    fopt
    [~,c,~] = objcon(xopt);
    c
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end

    % ------------Stalk geometry functions---------------
    function [r] = rpts(N,theta,dmaj,dmin,asymmetry,notch)
        r = zeros(1,N);
        for i = 1:N
            r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
                + ((dmaj/2)*sin(theta(i)))^2) + asymmetry(i) - notch(i);
        end
    end

    function [notch] = notch_fn(N,ndepth,nwidth,nloc,theta)
        notch = zeros(1,N);
        for i = 1:N
            notch(i) = ndepth/cosh((10/nwidth)*(theta(i)-nloc))^2;
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
