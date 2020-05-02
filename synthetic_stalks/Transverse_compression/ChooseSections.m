function ChooseSections(method,stalknums,dist,Table,error_indices,npoints,SaveName)
% FILENAME: ChooseSections.m
% AUTHOR: Ryan Larson
% DATE: 7/3/2019
%
% PURPOSE: The slice locations were created by the CT scan process, so they
% aren't all at "nice" distances. This function selects the desired
% cross-sections based on a method. 
% 
% 
% INPUTS:
%       method: String that determines how the cross-sections will be
%       selected. There are two methods: 'samedist' and 'all'. However, the
%       way later code was constructed means the only method that was used
%       was 'samedist'.
% 
%       stalknums: A row vector of integers corresponding to the stalk
%       numbers. These integers will be somewhere between 1 and 980.
% 
%       dist: A number (integer or float) that determines the desired slice
%       location, relative to the node. It can be positive or negative, and
%       is measured in mm. Maximum values should be approximately +/- 40mm.
% 
%       Table: A table containing the downsampled, centered, and rotated
%       cross-section boundaries. This is the table that is output from
%       Jared's script. Usually this will be Stalk_TableDCR, which loads as
%       a variable in StalksDCR_360pts_V2.mat.
% 
%       error_indices: Indices that had problems in pre-processing and
%       should be ignored. Also loads as a variable in
%       StalksDCR_360pts_V2.mat.
% 
%       npoints: The number of sample points the cross-sections have been
%       downsampled to. This was previously determined by Jared's
%       preprocessing script, and is only here to carry over the proper
%       indexing.
% 
%       SaveName: String name to save the data to as a .mat file. Must
%       include .mat in the name.
%       
% OUTPUTS:
%       N/A
%
%
% NOTES: 
%       Other outputs:
%           - Named .mat file that contains all the ellipse fit and PCA
%           data, which will be used when running transverse_wrapper_V4.m.
%       
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

allrows = size(Table,1);

switch method
    
    % Choose a number of cross-sections that are all at the same distance
    % from the node
    case 'samedist'
        indices = zeros(1,length(stalknums));

        
        % Step through stalks of interest (defined by range values)
        stalk = 1;
        for i = stalknums
            A = (Table.StkNum == i);
            
            % Get index of first row that is part of the current stalk
            for j = 1:length(A)
                if A(1) == 1
                    idx_first = 1;
                    break
                elseif A(j) == 1 && A(j-1) == 0
                    idx_first = j;
                end
            end
            
            % Get index of last row that is part of the current stalk
            for j = 2:length(A)
                if A(j) == 0 && A(j-1) == 1
                    idx_last = j-1;
                    break
                end
            end
            
            % Iterate through the indices that are part of the current
            % stalk and determine the index that is closest to the desired
            % value
            diffs = NaN(size(Table.StkNum));            
            for j = idx_first:idx_last
                diffs(j) = dist - Table.SlP(j);
            end
            
            % Get index of closest slice
            [~,closestIndex] = min(abs(diffs));
            
            % Add the index to indices output
            indices(stalk) = closestIndex;
            stalk = stalk + 1;
        end
        
        % Create table from indices for later reference
        selectedTable = Table(indices,:);
        
        % Get rid of any cross-sections that were chosen that are also
        % listed in error_indices
        for i = 1:length(error_indices)
            if ismember(error_indices(i),selectedTable.StkNum)
                row = find(selectedTable.StkNum == error_indices(i));
                selectedTable(row,:) = [];
            end
        end
        
        % Also get rid of any cross-sections that are made up of zeros
        % (some sort of detection error)
        deletesections = [];
        for i = 1:size(selectedTable,1)
            exterior_X = cell2mat(selectedTable.Ext_X(i));
            exterior_Y = cell2mat(selectedTable.Ext_Y(i));
            exterior_Rho = cell2mat(selectedTable.Ext_Rho(i));
            interior_X = cell2mat(selectedTable.Int_X(i));
            interior_Y = cell2mat(selectedTable.Int_Y(i));
            interior_Rho = cell2mat(selectedTable.Int_Rho(i));

            if (all(exterior_X == 0) || all(exterior_Y == 0) ||...
                    all(exterior_Rho == 0) || all(interior_X == 0) ||...
                    all(interior_Y == 0) || all(interior_Rho == 0))
                    
                deletesections = [deletesections, i];    
                msg = sprintf('Section %d deleted',i);
                disp(msg);
            end
        end
        
        %GET RID OF THE BAD SECTIONS HERE BY INDEXING
        selectedTable(deletesections,:) = [];
        
        % Save compiled slices in arrays for downstream use
        ext_X =     makearray(selectedTable,'Ext_X',npoints);
        ext_Y =     makearray(selectedTable,'Ext_Y',npoints);
        int_X =     makearray(selectedTable,'Int_X',npoints);
        int_Y =     makearray(selectedTable,'Int_Y',npoints);
        ext_T =     makearray(selectedTable,'Ext_T',npoints);
        ext_Rho =   makearray(selectedTable,'Ext_Rho',npoints);
        int_T =     makearray(selectedTable,'Int_T',npoints);
        int_Rho =   makearray(selectedTable,'Int_Rho',npoints);
        avg_rind_thick = selectedTable.rind_t;
        
        % Output all variables into mat file
        FolderName = pwd;
        SaveFile = fullfile(FolderName, SaveName);
        save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_T','ext_Rho',...
            'int_T','int_Rho','avg_rind_thick','indices','selectedTable',...
            'npoints','deletesections');

        
    % Choose all cross-sections (not really useful in hindsight, but was
    % written before full PCA process was defined)
    case 'all'        
        % Save compiled slices in arrays for downstream use
        ext_X =     makearray(Table,'Ext_X',npoints);
        ext_Y =     makearray(Table,'Ext_Y',npoints);
        int_X =     makearray(Table,'Int_X',npoints);
        int_Y =     makearray(Table,'Int_Y',npoints);
        ext_T =     makearray(Table,'Ext_T',npoints);
        ext_Rho =   makearray(Table,'Ext_Rho',npoints);
        int_T =     makearray(Table,'Int_T',npoints);
        int_Rho =   makearray(Table,'Int_Rho',npoints);
        avg_rind_thick = Table.rind_t;
        
        % Output all variables into mat file
        FolderName = pwd;
        SaveFile = fullfile(FolderName, SaveName);
        save(SaveFile,'ext_X','ext_Y','int_X','int_Y','ext_T','ext_Rho',...
            'int_T','int_Rho','avg_rind_thick','npoints');
        
    otherwise
        disp('Unknown method.');
end


end