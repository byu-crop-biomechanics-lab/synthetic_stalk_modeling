function [exported_array] = makearray(Table,Variable,npoints)
% FILENAME: makearray.m
% AUTHOR: Ryan Larson
% DATE: 7/3/2019
%
% PURPOSE: Convenience function found in some higher-level functions.
% Converts Table data for a given variable into an array, which is what
% downstream functions expect.
% 
% 
% INPUTS:
%       Table: Original data table (variable name)
% 
%       Variable: Variable name (string)
% 
%       npoints: Number of sample points around the cross-section
%       
% OUTPUTS:
%       exported_array: Array version of variable data in Table 
%
% NOTES:
%       
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

    num_slices = size(Table,1);
    exported_array = NaN(npoints,num_slices);
    
    if iscell(Table.(Variable))
        for i = 1:num_slices
            exported_array(:,i) = cell2mat(Table.(Variable)(i));
        end
    else
        for i = 1:num_slices
            exported_array(i,:) = Table.(Variable)(i);
        end
    end
        

end
