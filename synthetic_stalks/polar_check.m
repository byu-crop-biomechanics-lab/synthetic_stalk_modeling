dmaj = 20;
dmin = 15;

a = dmaj/2;
b = dmin/2;

N = 100;
theta = linspace(0,2*pi,N);

ndepth = 4;
nwidth = 5;
nloc = 3;

r = zeros(1,N);

for i = 1:N
    notch = ndepth/(cosh((10/nwidth)*(theta(i) - nloc))^2);
    r(i) = (a*b)/sqrt((b*cos(theta(i)))^2 + (a*sin(theta(i)))^2) - notch;
end


polarplot(theta,r)