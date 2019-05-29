%% getpolarfits.m:
% RL - 4/11/2019
% Assess how often poorly-fit outliers occur with the polar optimization
% routine

N = 100;
fits = zeros(N,1);

for i = 1:N
    [xopt, fopt, exitflag, output] = stalk_cross_fit_polar(i);
    fits(i) = fopt;
end

histogram(fits)

threshold = 0.06;

bad = sum(fits > threshold)
findbad = find(fits>threshold);
good = sum(fits <= threshold)

percentgood = good/N