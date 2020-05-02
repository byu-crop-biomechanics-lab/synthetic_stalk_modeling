function GenerateAllSTLs(Stalk_TableDCR,SESTable,npts,nTPCs,nalphaPCs)
% FILENAME: lookupStalk.m
% AUTHOR: Ryan Larson
% DATE: 2/7/2020
%
% PURPOSE: 
%   Wrap lookupStalk.m and GeneratePCAStalkSTL.m into a single function
%   that generates all the required STL files and names them accordingly
% 
% 
% INPUTS:
%       - Table: Table from Jared's SMALL_1500...something .mat
%       file
%       - set:
%       - hybrid: 
%       - density:
%       - replicate: 
%       - location: 
%       - stalk: Stalk number within the specified set
%       
%       
% OUTPUTS:
%       - stalknum: Stalk number in terms of the ~980 used for PCA
%
%
% NOTES: 
%       - 
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

Nstalks = size(SESTable,1);

for i = 1:Nstalks
    
    set = SESTable.sets(i);
    hybrid = SESTable.hybrids(i);
    density = SESTable.densities(i);
    replicate = SESTable.replicates(i);
    location = SESTable.locations(i);
    stalk = SESTable.stalks(i);
    
    stalknum = lookupStalk(Stalk_TableDCR,set,hybrid,density,replicate,location,stalk);
    
    for j = 1:nTPCs
        
        for k = 1:nalphaPCs
            
            FileName = sprintf('Stalk%d_%d_%d_PC_%d_%d',stalknum,set,stalk,j,k); % Change this to remove .stl, and then append .stl with Rind and Pith to make two distinct models
            
            GeneratePCAStalkSTL(stalknum,npts,nTPCs,nalphaPCs,FileName);
            
        end
    end

end

end