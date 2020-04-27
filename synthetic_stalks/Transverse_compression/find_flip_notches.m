function find_flip_notches(ChooseSectionsOutput,SaveName)
% FILENAME: find_flip_notches.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: Manually identify the cross-sections that need to be flipped 180
% degrees so the notch is on the left side. Cycles through the data
% contained in a .mat file output from ChooseSections.m, and allows the
% user to tag cross-sections that need to be flipped by entering "1" and
% then Enter. All other sections can be skipped through by hitting Enter.
% 
% 
% INPUTS:
%       ChooseSectionsOutput - A .mat file produced by ChooseSections.m
%       
% OUTPUTS:
%       flip_sections - A vector of 1s and empties (or non-1s) that must be
%       fed into flip_notches.m.
%
%
% NOTES: The reason find_flip_notches.m and flip_notches.m are not combined
% into a single function is so the user has the opportunity to correct for
% miskeyed cross-sections. flip_sections can be manually edited before
% being fed into flip_notches in case the user double-typed or had some
% other problem, since find_flip_notches.m only marches forward through the
% cross-sections, with no opportunity to make corrections during the
% process.
% 
% 
% VERSION HISTORY:
% V1 - Made into a function that works with the updated process flow
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load(ChooseSectionsOutput);

R_ext = ext_Rho;
R_int = int_Rho;
T = ext_T(:,1)';

N = size(R_ext,2);
flip_sections = zeros(N,1);

for i = 1:N
    % Plot each cross section to see if it needs to be flipped 180 degrees
    polarplot(T,R_ext(:,i));
    i
    s = input('Enter 1 if cross section needs to flip: ');
    if s ~= 1
        s = 0;
    end
    flip_sections(i) = s;
%     pause(); 
end
close;

% Save data as mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, SaveName);
save(SaveFile,'flip_sections');


end