function transverse_wrapper()
% FILENAME: transverse_wrapper.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: Wrap the majority of the data production process into a single
% script
% 
% 
% INPUTS:
%      
%       
% OUTPUTS:
%       
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
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

%% Process
% Load in the 360-point DCR version of Jared's table
load StalksDCR_360pts.mat

% Choose the cross-section samples
SliceOutputName = 
ChooseSections('samedist',range,Stalk_TableDCR,npoints,SliceOutputName)

% Manually find the cross-sections that need to be flipped 180 degrees
FlipName
find_flip_notches(ChooseSectionsOutput,FlipName)

% Flip the cross-sections that need to be flipped, according to the vector
% of flip indicators
flip_notches(Flip_Indices,ChooseSectionsOutput,FlippedOutputName)

% Calculate the ellipse fits
ellipse_fitting_V2(FileName,SaveName)

% Plot the ellipse fits and see if any of them have problems
for i = 1:10
    
    
end


good = input('Are all the ellipse fits good?');
% Might need an automatic way to list the slices that are bad

if good == 1
    % Run PCA
    PCA_ellipse_fits(FileName,SaveName)
else
    % Remove the problem ellipses and then run PCA again
    remove_problem_ellipses(OriginalEllipseFits,problem_indices,GoodEllipseFits)
    
    % Run PCA
    PCA_ellipse_fits(FileName,SaveName)
end


% If the PCA is now good, create the Abaqus Python scripts
create_cases(NEPCdata,EllipseData,numNEPCs)


end