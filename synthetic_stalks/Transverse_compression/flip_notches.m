function flip_notches(Flip_Indices,ChooseSectionsOutput,FlippedOutputName)
% FILENAME: find_flip_notches.m
% AUTHOR: Ryan Larson
% DATE: 6/18/19
%
% PURPOSE: Take the output of find_flip_notches.m and rotate cross-sections
% that were manually tagged. This corrects errors in the automatic
% preprocessing that aren't caught by Jared's script.
% 
% 
% INPUTS:
%       Flip_Indices: The .mat file that contains flip_sections, which is a
%       vector of 1s and 0s that identifies the cross-sections needing to
%       be flipped. 
%
%       ChooseSectionsOutput - A .mat file produced by ChooseSections.m
%
%       FlippedOutputName - The name for the output .mat file. Make sure to
%       end the name with FLIPPED for consistency.
%       
% OUTPUTS:
%       N/A
%
% NOTES:
%       
% 
% 
% VERSION HISTORY:
% V1 - Made into a function that works with the updated process flow
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

load(ChooseSectionsOutput);
load(Flip_Indices);

ext_rho = [];
int_rho = [];
T = ext_T(:,1)';
N = size(ext_X,2);

flipped = []; % Hold the indices of the flipped sections
for i = 1:N
    if flip_sections(i) == 1 
        flipped = [flipped; i];
        
        % Rotate and reorder external and internal points
        [~, ~, ~, ~, ~, ~, xp_ext, yp_ext, ~, ~] = reorder_V2(ext_X(:,i), ext_Y(:,i), pi);
        [~, ~, ~, ~, ~, ~, xp_int, yp_int, ~, ~] = reorder_V2_interior(int_X(:,i), int_Y(:,i), pi, 0, 0);
        
        [~, ~, x_ext, y_ext, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_ext, yp_ext, 0);
        [~, ~, x_int, y_int, ~, ~, ~, ~, ~, ~] = reorder_V2_interior(xp_int, yp_int, 0, 0, 0);
        
        % Redefine the appropriate row in the main XY data
        ext_X(:,i) = x_ext;
        ext_Y(:,i) = y_ext;
        int_X(:,i) = x_int;
        int_Y(:,i) = y_int;
    end
end

% Convert data from Cartesian to polar
size(T)
for i = 1:N
    for j = 1:length(T)
        ext_Rho(j,i) = sqrt(ext_X(j,i)^2 + ext_Y(j,i)^2);
        int_Rho(j,i) = sqrt(int_X(j,i)^2 + int_Y(j,i)^2);
    end
end

% Copy the variable selectedTable as a basis for the new data
flippedTable = selectedTable;

% Apply orientation corrections to the appropriate rows in flippedTable
for i = 1:N
    flippedTable.Ext_X{i} = ext_X(:,i);
    flippedTable.Ext_Y{i} = ext_Y(:,i);
    flippedTable.Ext_Rho{i} = ext_Rho(:,i);
    flippedTable.Int_X{i} = int_X(:,i);
    flippedTable.Int_Y{i} = int_Y(:,i);
    flippedTable.Int_Rho{i} = int_Rho(:,i);
    
end

% Save data as mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, FlippedOutputName);
save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_Rho','int_Rho','ext_T',...
    'int_T','avg_rind_thick','flip_sections','flippedTable','npoints');

end