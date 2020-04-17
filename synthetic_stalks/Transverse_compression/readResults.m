function [ResultsCell] = readResults(ResultsTxtFile,ResultsSaveFile)

% This is designed for reading results of results files that come from
% transverse_wrapper_V4.m (not the sensitivity study).

ResultsCell = {};

s = tdfread(ResultsTxtFile,'\t');

Group = s.Group;
ID = s.Group;
Case = s.Case;
RFx = s.RFx;
RFy = s.RFy;

groups = transpose(unique(Group));

for group = groups
    group
    inds = Group(:,1) == group;
    Results = [ID(inds),Case(inds),RFx(inds),RFy(inds)];
    
    ResultsCell(1,group) = {Results};    
end

% Output all variables into mat file
FolderName = pwd;
SaveFile = fullfile(FolderName, ResultsSaveFile);
save(SaveFile,'ResultsCell');

end