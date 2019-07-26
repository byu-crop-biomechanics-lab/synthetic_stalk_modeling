
N = 2000;

% Equations
deflection = 0.5;
E = 10000;
V = 100;

D = zeros(N,1);
d = zeros(N,1);
L = zeros(N,1);
stress = zeros(N,1);

% Generate variable values that meet conditions
i = 1;
iteration = 0;

while i <= N
    iteration = iteration + 1;
    D(i) = unifrnd(0.1,25);
    d(i) = unifrnd(0.1,25);
    
    t = (D(i)-d(i))/2;
    Area = (pi/4)*(D(i)^2 - d(i)^2);
    L(i) = V/Area;
    
    if t < 0.1 || d(i) >= D(i) || L(i) < 5 || L(i) > 75
        continue
    end
    
    % Calculate max bending stress
    I = (pi/64)*(D(i)^4 - d(i)^4);
    F = (48*deflection*E*I)/L(i)^3;
    M = F*L(i)/2;
    stress(i) = M*D(i)/(2*I);
    
    i = i + 1;
    
end

noise = unifrnd(-100,100,size(stress));

outputs = [stress+noise, D, d, L];

figure(1);
scatter3(D,d,L,'.');
xlabel('D');
ylabel('d');
zlabel('L');
title('Variable trends');

figure(2);
scatter3(D,d,outputs(:,1),'.');
xlabel('D');
ylabel('d');
zlabel('stress');
title('Stress response');

