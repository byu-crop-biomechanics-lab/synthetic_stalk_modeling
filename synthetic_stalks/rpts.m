function [r] = rpts(N,theta,dmaj,dmin,asymmetry,notch)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2) + asymmetry(i) - notch(i);
    end
end