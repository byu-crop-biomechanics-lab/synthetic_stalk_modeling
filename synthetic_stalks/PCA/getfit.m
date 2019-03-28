function [fitmetric] = getfit(xreal,xsynth,yreal,ysynth,zreal,zsynth)
    % Get the overall fit of a synthetic shape against a real shape, in a
    % single value similar to a standard deviation
    
    N = length(xreal);
    
    if nargin == 2
        yreal = zeros(N);
        ysynth = zeros(N);
        zreal = zeros(N);
        zsynth = zeros(N);
    
    elseif nargin == 4
        zreal = zeros(N);
        zsynth = zeros(N);
        
    elseif mod(nargin,2) == 1
        error('Odd number of arguments');
    end
    
    sum = 0;
    
    for i = 1:N
        xsq = (xsynth(i) - xreal(i))^2;
        ysq = (ysynth(i) - yreal(i))^2;
        zsq = (zsynth(i) - zreal(i))^2;
        sum = sum + xsq + ysq + zsq;
    end
    
    fitmetric = sqrt(sum/N);
    
end