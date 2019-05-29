function [ ext_X, ext_Y, xbar, ybar] = exterior_boundaries_V4( I, thresh, plotting )

%====================================================================================
% FILE:   exterior_boundaries.m
% AUTHOR: Douglas Cook 
% 
% PURPOSE: This function identifies the exterior boundaries of a CT scan of a plant cross-section
%
%
% INPUTS: 
%         I       - a 16-bit grayscale image containing a stalk cross-section.  
%         thresh  - the grayscale threshhold for identifying the stalk boundaries (between 0 and 65535 for 16-bit grayscale images)
% 
% OUTPUTS: 
%       ext_X and ext_Y - the x and y coordinates of the exterior boundary between the stalk and the air. 
%       xbar, ybar      - the x and y coordinates of the geometric center of the boundary curve
% 
%
% VERSION HISTORY
%       Version R2: rind_thickness was removed as an output since a new function "avg_rind_thickness_normal_method()" 
%       now produces a much more reliable estimation of rind thickness.
% 
%       R3 - all int_X, int_Y outputs removed along with associated code since it is no
%       longer used.
%
%       V4 - 3/17/18 - code cleaned up and unnecessary/obsolete portions removed.
%                    - plotting option added
%
% NOTES: 
%       
%


% CONTOUR ANALYSIS
[C] = contourc(double(I), [1 1]*thresh);    % find the contours at threshhold level
index = C(1,:) == thresh;                   % find "thresh" in the first row of the contour matrix - these correspond to gaps between curves
index = C(2,index);                         % get the lengths of each curve in the contour matrix
index = sort(index, 'descend');             % sort the curves by length


k = 1;          
while index(k) > 400                                    % loop through all contours that are longer than 500 pixels long
    
    if index(k) == index(k+1)                           % in the rare event that two curves have the same length, use find(_____,1) to get the first, find(____,2) to get the second.   
        c_st(k) = find(C(2,:)==index(k),1) + 1;         % the starting location of the curve in the contour matrix
        c_end(k) = c_st(k) + index(k)-1;                % the ending location 
        curve{k} = C(:,c_st(k):c_end(k));               % a cell array containing each contour curve data 
        diam(k) = max(max(curve{k}')-min(curve{k}'));   % calculates the maximum x or y span of the curve
        k  = k + 1;
        
        temp = find(C(2,:)==index(k),2) + 1;            % the starting location of the interior curve in the contour matrix
        c_st(k) = temp(2);
        c_end(k) = c_st(k) + index(k)-1;                % the ending location 
        curve{k} = C(:,c_st(k):c_end(k));               % each contour curve data 
        diam(k) = max(max(curve{k}')-min(curve{k}'));   % calculates the maximum x or y span of the curve
        k  = k + 1;
    else
        c_st(k) = find(C(2,:)==index(k)) + 1;           % the starting location of the curve in the contour matrix
        c_end(k) = c_st(k) + index(k)-1;                % the ending location 
        curve{k} = C(:,c_st(k):c_end(k));               % each contour curve data 
        diam(k) = max(max(curve{k}') - min(curve{k}')); % calculates the maximum x or y span of the curve
        k = k + 1;
    end
end


[~,dindex] = sort(diam,'descend');              % sort curves by diameter

exterior = curve{dindex(1)};                    % the exterior has the largest diameter
%interior = curve{dindex(2)};                   % the interior curve is no longer used.

% external and internal contours as column vectors 
ext_X = exterior(1,:)';                         % take the first row and tanspose
ext_Y = exterior(2,:)';                         % second row and transpose
%int_Xi = interior(1,:)';
%int_Yi = interior(2,:)';

% comparisons show that the mean of each coordinate is within 1 pixel of the true geometric center.
xbar = mean(ext_X);
ybar = mean(ext_Y);



% PLOTTTING (OPTIONAL)

if plotting == 1 
    figure;
    imshow(I)
    hold on
    plot(ext_X,ext_Y,'y')
end

end