%% Setup
% Two varieties each of wheat and barley, and four cultivars of corn stover
% were used.

% For wheat and barley, the second internode from the top of the plant was
% examined.

% For corn, the internodes above and below the cob location were examined.
% 

%% Bending and Tension Tests
% Wheat and barley had 1 specimen per internode (1 internode examined per
% stalk)
% Corn was not examined

%% Compression Tests:
% Wheat and barley: 2 samples cut per internode (1 internode examined per
% stalk)
% "7 corn specimens for each variety were cut from the same internodal
% region, one internode above and below the cob location"
% I think this means the 7 specimens were spread across the two internodes
% of interest for corn?

%% Setting up data
clc
clear all
close all

% Given data for means (all in GPa, and the vector order is the same as in
% Table 3):
bend3_mean = [2.2; 1.1; 1.3; 1.1];
bend4_mean = [2.2; 1.3; 1.4; 1.3];
compression_mean = [0.60; 0.90; 0.42; 0.97; 0.26; 0.38; 0.47; 0.40];
tension_mean = [7.3; 4.9; 3.4; 3.5];

% Given data for standard error of means (all in GPa, with Table 3 order):
bend3_sem = [0.23; 0.076; 0.087; 0.063];
bend4_sem = [0.18; 0.020; 0.058; 0.035];
compression_sem = [0.11; 0.47; 0.08; 0.38; 0.06; 0.04; 0.12; 0.10];
tension_sem = [0.92; 0.74; 0.57; 0.37];

% % Guessing sample sizes to calculate standard deviations
% wstalks = 200;   % # of wheat stalks sampled
% bstalks = 20;   % # of barley stalks sampled
% cstalks = 20;   % # of corn stalks sampled

% 7 samples were taken per variety for any given test
n_w_c = 7;
n_w_bt = 7;
n_b_c = 7;
n_b_bt = 7;
n_c_c = 7;

%% Calculate the synthetic standard deviation for each test:

bend3_stdev = zeros(length(bend3_sem),1);
bend4_stdev = zeros(length(bend4_sem),1);
compression_stdev = zeros(length(compression_sem),1);
tension_stdev = zeros(length(tension_sem),1);

for i = 1:length(bend3_sem)
    % Wheat cases, all except compression
    if i <= 2
        bend3_stdev(i) = bend3_sem(i)*sqrt(n_w_bt);
        bend4_stdev(i) = bend4_sem(i)*sqrt(n_w_bt);
        tension_stdev(i) = tension_sem(i)*sqrt(n_w_bt);
    % Barley cases, all except compression
    elseif i <= 4
        bend3_stdev(i) = bend3_sem(i)*sqrt(n_b_bt);
        bend4_stdev(i) = bend4_sem(i)*sqrt(n_b_bt);
        tension_stdev(i) = tension_sem(i)*sqrt(n_b_bt);
    end
end

for i = 1:length(compression_sem)
    % Wheat cases, compression
    if i <= 2
        compression_stdev(i) = compression_sem(i)*sqrt(n_w_c);
    % Barley cases, compression
    elseif i <= 4
        compression_stdev(i) = compression_sem(i)*sqrt(n_b_c);
    % Corn cases (compression only)
    else
        compression_stdev(i) = compression_sem(i)*sqrt(n_c_c);
    end
end

%% Create synthetic data points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now that we have the synthetic standard deviations, we can create
% synthetic data as follows:

testlabels = {'3-Point Bending','4-Point Bending','Compression','Tension'};

%% Amidon wheat
aw3 = bend3_stdev(1)*randn(n_w_bt,1) + bend3_mean(1);
% % sanity check:
% stats = [mean(aw3) std(aw3)]
aw4 = bend4_stdev(1)*randn(n_w_bt,1) + bend4_mean(1);
awc = compression_stdev(1)*randn(n_w_c,1) + compression_mean(1);
awt = tension_stdev(1)*randn(n_w_bt,1) + tension_mean(1);

Amidonwheat = [aw3, aw4, awc, awt];
figure(1);
boxplot(Amidonwheat,'Labels',testlabels);
title('Amidon (wheat)');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');

%% Westbred 936 wheat
wb3 = bend3_stdev(2)*randn(n_w_bt,1) + bend3_mean(2);
wb4 = bend4_stdev(2)*randn(n_w_bt,1) + bend4_mean(2);
wbc = compression_stdev(2)*randn(n_w_c,1) + compression_mean(2);
wbt = tension_stdev(2)*randn(n_w_bt,1) + tension_mean(2);

Westbred936 = [wb3, wb4, wbc, wbt];
figure(2);
boxplot(Westbred936,'Labels',testlabels);
title('Westbred 936 (wheat)');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');

%% Bowman barley
bb3 = bend3_stdev(3)*randn(n_w_bt,1) + bend3_mean(3);
bb4 = bend4_stdev(3)*randn(n_w_bt,1) + bend4_mean(3);
bbc = compression_stdev(3)*randn(n_w_c,1) + compression_mean(3);
bbt = tension_stdev(3)*randn(n_w_bt,1) + tension_mean(3);

Bowman = [bb3, bb4, bbc, bbt];
figure(3);
boxplot(Bowman,'Labels',testlabels);
title('Bowman (barley)');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');

%% Fragile Stem 1 barley
fs3 = bend3_stdev(4)*randn(n_w_bt,1) + bend3_mean(4);
fs4 = bend4_stdev(4)*randn(n_w_bt,1) + bend4_mean(4);
fsc = compression_stdev(4)*randn(n_w_c,1) + compression_mean(4);
fst = tension_stdev(4)*randn(n_w_bt,1) + tension_mean(4);

FragileStem = [fs3, fs4, fsc, fst];
figure(4);
boxplot(FragileStem,'Labels',testlabels);
title('Fragile Stem 1 (barley)');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');

%% Create bar plot with error bars %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Amidon wheat
awmeans = [bend3_mean(1),bend4_mean(1),compression_mean(1),tension_mean(1)];
awerrors = [bend3_sem(1),bend4_sem(1),compression_sem(1),tension_sem(1)];
figure(5);
barwitherr(awerrors,awmeans);
set(gca,'XTickLabel',{'3-Point Bending','4-Point Bending','Compression','Tension'})
title('Amidon wheat - Means and Standard Error');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');

% Westbred 936 wheat
wbmeans = [bend3_mean(2),bend4_mean(2),compression_mean(2),tension_mean(2)];
wberrors = [bend3_sem(2),bend4_sem(2),compression_sem(2),tension_sem(2)];
figure(6);
barwitherr(wberrors,wbmeans);
set(gca,'XTickLabel',{'3-Point Bending','4-Point Bending','Compression','Tension'})
title('Westbred 936 wheat - Means and Standard Error');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');

% Bowman barley
bbmeans = [bend3_mean(3),bend4_mean(3),compression_mean(3),tension_mean(3)];
bberrors = [bend3_sem(3),bend4_sem(3),compression_sem(3),tension_sem(3)];
figure(7);
barwitherr(bberrors,bbmeans);
set(gca,'XTickLabel',{'3-Point Bending','4-Point Bending','Compression','Tension'})
title('Bowman barley - Means and Standard Error');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');

% Fragile Stem 1 barley
fsmeans = [bend3_mean(4),bend4_mean(4),compression_mean(4),tension_mean(4)];
fserrors = [bend3_sem(4),bend4_sem(4),compression_sem(4),tension_sem(4)];
figure(8);
barwitherr(bberrors,bbmeans);
set(gca,'XTickLabel',{'3-Point Bending','4-Point Bending','Compression','Tension'})
title('Fragile Stem 1 barley - Means and Standard Error');
xlabel('Testing Method');
ylabel('Modulus of Elasticity (GPa)');
