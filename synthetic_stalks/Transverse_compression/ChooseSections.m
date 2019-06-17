function ChooseSections(method,range,Table,SaveName)
% ChooseSections.m: Determine the cross-sections to compile, which is
% determined by a method

% range: For Stalk_Table, this must be a row vector of two integer values
% from 1 to 990.

allrows = size(Table,1);

switch method
    % Choose a number of cross-sections that are all at the same distance
    % from the node
    case 'samedist'
        indices = zeros(1,(range(2) - range(1) + 1));
        dist = input('Choose approximate slice distance to use: ');
        
        % Step through stalks of interest (defined by range values)
        stalk = 1;
        for i = range(1):range(2)
            A = (Table.StkNum == i); % NEED AN INTELLIGENT WAY TO GET THE STARTING AND ENDING INDICES FOR THE PARTS OF A THAT ARE TRUE (1)
            
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
        
        % Save compiled slices in arrays for downstream use
        ext_X =     makearray(selectedTable,'Ext_X',360);
        ext_Y =     makearray(selectedTable,'Ext_Y',360);
        int_X =     makearray(selectedTable,'Int_X',360);
        int_Y =     makearray(selectedTable,'Int_Y',360);
        ext_T =     makearray(selectedTable,'Ext_T',360);
        ext_Rho =   makearray(selectedTable,'Ext_Rho',360);
        int_T =     makearray(selectedTable,'Int_T',360);
        int_Rho =   makearray(selectedTable,'Int_Rho',360);
        avg_rind_thick = Table.rind_t(range(1):range(2));
        
        
        % Output all variables into mat file
        FolderName = pwd;
        SaveFile = fullfile(FolderName, SaveName);
        save(SaveFile);
        
    case 'wholestalk'
        % Choose a range of stalk numbers, and all the slices from each of
        % the chosen stalks will be chosen
        
    otherwise
        disp('Unknown method.');
end


end


function [exported_array] = makearray(Table,Variable,npoints)
    % Make an array of the chosen variable in the format required to work
    % with the downstream process in ellipse_fittingV1.m and
    % PCA_ellipse_fits.m
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