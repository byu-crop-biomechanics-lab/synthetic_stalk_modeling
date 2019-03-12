% FILENAME: notch_finder
% AUTHOR: Aaron Lewis
% DATE: 8/17/18
%
%
% PURPOSE: Identify the location of the notch on the cross-section of a
%          stalk of corn. 
% 
% 
% INPUTS: 
% 
% 
% OUTPUTS:        
%
%
% NOTES: This is done by identifying the peaks on a polar curve
%        representation of the exterior boundary of the cross-section. The two
%        peaks that are closest together will be the upper and lower boundaries of
%        the notch location. From there, the notch is simply the minimum value of
%        that range.
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
% 
% -------------------------------------------------------------------------

clc
clear


% Inputs
% filename = '20180227_2InterNodal Y_S1.mat';
filename = 'RealData_DataSet_2.mat';
npoints = 100;

[s, xi, yi, ti, rhoi] = rotate_and_reorder( filename, npoints );        % Rotates, reorders, and downsamples the stalk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npoints = length(ti);

newrho = [rhoi rhoi];                                            % Duplicates the rho vector, just in case the notch is at the very start or very end of the rho curve
sel = (max(newrho)-min(newrho))/32;                              % I played around with this number until it identified all of the peaks of the rho curve.
                                                                 %     The smaller the number,
                                                                 %     the more sensitive it would
                                                                 %     be to the smaller peaks
peakfinder(newrho,sel);                                        % Graphs the peakfinder (for convenience)
[peaklocs, peakmags] = peakfinder(newrho,sel);                               % Grabs the indecies of each peak on the rho curve

peakmags
[~, highlocs] = sort(peakmags,'descend')

cutoff = (max(newrho) + min(newrho))/2;

for m = 1:length(peakmags)
    if peakmags(m) < cutoff
        peakmags(m) = NaN;
    end
end
I = peakmags>cutoff
peakmags
    



% The way this function identifies the notch location is by finding the two peaks of
% the rho curve that are closest together. This seems to reliably give an
% upper and lower bound around where the center of the notch is.
for i = 2:length(peaklocs);
    peakdifs(i-1) = peaklocs(i) - peaklocs(i-1);                 % Finding the distances between each pair of peaks
end

[~, index] = min(peakdifs);                                      % Identifies the index of the smallest peak-to-peak distance
notch_range = newrho(peaklocs(index):peaklocs(index+1));         % Creates a vector from the two closest peaks (the notch will be the minimum of this vector)
[notch, notchloc] = min(notch_range);                            % Finds the index of the notch
notchloc = notchloc-1;                                           % The index will be 1 too many
notchloc = peaklocs(index) + notchloc;                           % Add the index of the notch to the first peak to get the complete index of the notch

% Since the rho vector was doubled up at the beginning of this process, its
% possible that the notch could be located in the second half. This can be
% fixed by simply subtracting the length of the rho vector from the index.
if notchloc > 100         
    notchloc = notchloc - npoints;
end


% % Plotting the evidence
% subplot(1,2,1)
% plot(ti,rhoi,ti(notchloc),rhoi(notchloc),'ro')
% subplot(1,2,2)
% plot(xi, yi, xi(notchloc), yi(notchloc), 'ro')
% axis equal











