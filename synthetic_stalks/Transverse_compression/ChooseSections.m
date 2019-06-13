function [indices] = ChooseSections(method,range)
% ChooseSections.m: Determine the cross-sections to compile, which is
% determined by a method

switch method
    % Choose a number of cross-sections that are all at the same distance
    % from the node
    case 'samedist'
        indices = zeros(1,(range(2)-range(1)));
        dist = prompt('Choose approximate slice distance to use: ');
        
            
            
            
            
        end
    otherwise
        disp('Unknown method.');
end


end