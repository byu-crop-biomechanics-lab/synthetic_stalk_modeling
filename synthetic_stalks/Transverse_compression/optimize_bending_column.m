function [xopt, fopt, exitflag, output] = optimize_bending_column()

    % ------------Starting point and bounds------------
    %vars D    d    E
    x0 = [19   14   10000];
    ub = [20   20   100000];
    lb = [0.1  0.1  100];

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        % Design variables
        D = x(1);
        d = x(2);
        E = x(3);
        
        % Other analysis variables
        deflection = 0.5;
        L = 20;
        
        % Analysis functions
        t = (D-d)/2;
        Area = (pi/4)*(D^2 - d^2);
        V = Area*L;
        I = (pi/64)*(D^4 - d^4);
        F = (48*deflection*E*I)/L^3;
        M = F*L/2;
        stress = M*D/(2*I);
        
        % Objective function
        f = stress; % Make objective negative to maximize value
        
        % Inequality constraints
        c = zeros(2,1);
        c(1) = d - D;           % d <= D
        c(2) = 0.05 - t;        % 0.05 <= t
%         c(3) = 2*rindt - dmaj;  % Rind thickness can't be larger than the cross-section
%         c(4) = 2*rindt - dmin;  % Rind thickness can't be larger than the cross-section
%         c(5) = 75000 - biomass;
%         c(6) = biomass - 100000;
        
        % Equality constraints
        ceq = V - 400;
%         ceq = [];

        

    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'display', 'iter-detailed');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end
end