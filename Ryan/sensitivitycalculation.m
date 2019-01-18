clear all
E = 1.503550000000000e+04;
Ibase = 4.6759e+05;
EIbase = E*Ibase;
I = [4.6992e+05; 4.7226e+05; 4.7460e+05; 4.7694e+05; 4.9097e+05]
EI = E*I
factor = [1.005; 1.01; 1.015; 1.02; 1.05]
deltabase = 0.0433;
delta = [0.0433; 0.0432; 0.0432; 0.0432; 0.0431]

sensitivity = zeros(length(I),1);
for i = 1:length(sensitivity)
    num = (abs(deltabase-delta(i))/deltabase)*100;
    den = (abs(EIbase - EI(i))/EIbase)*100;
    sensitivity(i) = num/den;
end
sensitivity