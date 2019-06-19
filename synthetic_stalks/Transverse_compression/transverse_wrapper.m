function transverse_wrapper(range,slicedist)
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
hold off
close all;
set(0,'DefaultFigureWindowStyle','docked');

r1 = num2str(range(1));
r2 = num2str(range(2));
dist_int = num2str(round(abs(slicedist),0));
deci = abs(slicedist) - round(abs(slicedist),0);
dist_deci = num2str(deci);
dist_deci = erase(dist_deci,'0.');

if slicedist > 0
    if strcmp(dist_deci,'0')
        slicepos = strcat('_Above_',dist_int);
    else
        slicepos = strcat('_Above_',dist_int,'_',dist_deci);
    end
elseif slicedist < 0
    if strcmp(dist_deci,'0')
        slicepos = strcat('_Below_',dist_int);
    else
        slicepos = strcat('_Below_',dist_int,'_',dist_deci);
    end
else
    slicepos = strcat('_At_Node');
end

output_prefix = strcat('Stalks_',r1,'_',r2,slicepos);

% Choose the cross-section samples
ChooseSectionsName = strcat(output_prefix,'_Sampled.mat');
ChooseSections('samedist',range,slicedist,Stalk_TableDCR,npoints,ChooseSectionsName)

% Manually find the cross-sections that need to be flipped 180 degrees
FlipName = strcat(output_prefix,'_flip_sections.mat');
find_flip_notches(ChooseSectionsName,FlipName)

while 1
    fixes_needed = input('Does the flip vector need manual correction? Y/N ','s');
    switch fixes_needed
        case 'Y'
            load(FlipName);
            openvar('flip_sections');
            disp('Giving control to keyboard for manual editing of flip variable.');
            disp('Use dbcont command to exit keyboard editing mode.');
            keyboard;
            break
        case 'N'
            break
            
        otherwise
            disp('Not a recognized answer. Please try again.');
    end
end


% Flip the cross-sections that need to be flipped, according to the vector
% of flip indicators
% load(FlipName);
FlippedOutputName = strcat(output_prefix,'_FLIPPED.mat');
flip_notches(FlipName,ChooseSectionsName,FlippedOutputName);

% Calculate the ellipse fits
EllipseName = strcat(output_prefix,'_Ellipses.mat');
ellipse_fitting_V2(FlippedOutputName,EllipseName);

% Plot the ellipse fits and see if any of them have problems
load(EllipseName);

problem_indices = [];
for i = 1:(range(2) - range(1) + 1)
    i
    polarplot(ELLIPSE_T(i,:),ELLIPSE_R_ext(i,:));
    hold on
    polarplot(ELLIPSE_T(i,:),ELLIPSE_R_int(i,:));
    hold off
    s = input('Enter 1 if the ellipse fit has a problem: ');
    s
    if s == 1
        problem_indices = [problem_indices, i];
    else
        continue
    end    
end

while 1
    fixes_needed = input('Does the ellipse problems vector need manual correction? Y/N ','s');
    switch fixes_needed
        case 'Y'
            load(FlipName);
            openvar('problem_indices');
            disp('Giving control to keyboard for manual editing of flip variable.');
            disp('Use dbcont command to exit keyboard editing mode.');
            keyboard;
            break
        case 'N'
            break
            
        otherwise
            disp('Not a recognized answer. Please try again.');
    end
end

problem_indices

NEPCName = strcat(output_prefix,'_PCA.mat');

if isempty(problem_indices)
    % Run PCA
    PCA_ellipse_fits(EllipseName,NEPCName);
    
    % Create the Abaqus Python scripts
    create_cases(NEPCName,EllipseName,ChooseSectionsName,problem_indices,5);
else
    % Remove the problem ellipses and then run PCA again
    GoodEllipseFits = strcat(output_prefix,'_GoodEllipses.mat');
    remove_problem_ellipses(EllipseName,problem_indices,GoodEllipseFits);
    
    % Run PCA
    PCA_ellipse_fits(GoodEllipseFits,NEPCName);
    
    % Create the Abaqus Python scripts
    create_cases(NEPCName,GoodEllipseFits,ChooseSectionsName,problem_indices,5);
end

set(0,'DefaultFigureWindowStyle','normal');

end