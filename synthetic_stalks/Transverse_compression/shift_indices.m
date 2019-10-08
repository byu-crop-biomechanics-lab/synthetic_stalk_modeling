function [newgoodstalknums] = shift_indices(problem_indices,Nstalks,goodstalknums)
    % Deal with removal of problem_indices by shifting and renumbering the 
    
    % goodstalknums is used as the chosen set of stalks for analysis, minus
    % the bad ones that are also found in problem_indices
    
    original_indices = linspace(1,Nstalks,Nstalks)';
    original_indices(problem_indices) = [];     % Remove the problem indices from the list of indices
    
    good_indices = original_indices;    % Vector of indices with skipped values called out by problem_indices
    
    shifts = zeros(size(good_indices));
    
    % Calculate shifts and map them to the good_indices
    expected = 1;
    for i = 1:length(good_indices)
        if good_indices(i) ~= expected
            diff = good_indices(i) - expected;
            shifts(i:end) = shifts(i:end) + diff - 1;
        end
        expected = good_indices(i);
    end
    
    % Shift goodstalknums according to the mapped shift values
    newgoodstalknums = zeros(size(goodstalknums));
    for i = 1:length(goodstalknums)
        shiftval = shifts(find(good_indices==goodstalknums(i)));
        newgoodstalknums(i) = goodstalknums(i) - shiftval;
    end

end