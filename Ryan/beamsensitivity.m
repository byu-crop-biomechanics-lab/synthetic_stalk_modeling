%% Setup %%
clear all
Nx = 3;     % Number of sections on the left, or x, side of the anvil
Ny = 3;     % Number of sections on the right, or y, side of the anvil
lx = [100 100 100 100];     % Section lengths in mm
ry = [100 100 100];       % Section lengths in mm

nodes = zeros(1,length(lx) + length(ry) + 1);
internodes = linspace(1,length(nodes)-1,length(nodes)-1);

F = 10; % N
Evalues = [14626 14840 15187 13284 16610 15666]; % N/mm^2, taken from 6
                                                 % models used in SES2018
                                                 % spreadsheet
E = mean(Evalues);  % We will vary EI only by adjusting I

for i = 1:length(lx)
    iremainx = (i+1):length(lx);   % Remaining indices on x side
    totx = 0;
    totx = totx + lx(i);
    for j = 1:length(iremainx)
        totx = totx + lx(iremainx(j));
    end
    nodes(i) = totx;
end

toty = 0;
for i = 1:length(ry)
    toty = toty + ry(i);
    nodes(i + 1 + length(lx)) = toty;
end
nodes

I = zeros(1,length(internodes));
maj = 30;   % mm, major diameter of ellipse approximation for cross section
min = 20;   % mm, minor diameter of ellipse approximation for cross section

I(1) = (pi/4)*min*maj^3  % This will be the I value at the far left node
Ipercentchange = 0.05;  % 5% change in I between internodes                            
for i = 2:length(I)
    I(i) = I(i-1)*(1+Ipercentchange);
end

factor = 1.05
I(3) = factor*I(3)

Ianvil = (I(length(lx)+1) + I(length(lx)))/2

%% Calculations
leftterms = 0;
rightterms = 0;
for e = 1:length(internodes)
    x1 = nodes(e);
    x2 = nodes(e+1);
    if x2 > x1
        rightterms = rightterms + (1/I(e))*(x2^3 -x1^3);
    elseif x1 > x2
        leftterms = leftterms + (1/I(e))*(x1^3 - x2^3);
    end
end
delta = (F/(3*E))*(rightterms + leftterms)
