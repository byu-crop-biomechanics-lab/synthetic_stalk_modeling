function [Inertia_Centroid_Y,Iuu] = Center_Inertia(Int_x,Int_y,Ext_x,Ext_y,tol)
% FILENAME: Center_Inertia.m
% AUTHOR: Michael Ottesen
% DATE: 5/13/2020
%
% PURPOSE: Calculate intertia center in y direction starting from gemetric
%          center.
% 
% INPUTS:
%       Data points: Exterior and interior data points needed for a single 
%       cross section. They must be in cartesian coordinates.
%       
% OUTPUTS:       
%       Inertia_Centroid_Y: The position if the inertia
%
% NOTES: 
%       - 
% 
% PSEUDO-CODE:
%   Repeat last point in data at end of each data array.
%   Divide rind section into top and bottom halves.
%   Get initial inertias and geometric values.
%   Make first guess (guessing endpoints).
%   Initialize variables for the 'while' loop.
%   Begin loop through each cross section.
%       Begin while loop for while inertias are not equal.
%           Check to see if intertias are within the tolerance.
%           Use a bisection method to narrow down correct inertia center.
%               Make guess (half way between the two guessing endpoints).
%               Calcualte new inertia values for top and bottom using the
%               parallel axis theorem.
%               Depending on the side the true center is, adjust the
%               endpoints to use new guess.
%       end while loop
%   Close for loop
%   
%        
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

% Repeat the first point at the end of the arrays
int_x = [Int_x,Int_x(1)];
int_y = [Int_y,Int_y(1)];
ext_x = [Ext_x,Ext_x(1)];
ext_y = [Ext_y,Ext_y(1)];

% Determine if number of points is odd or even (m is a place holder for
% getting the top and bottom sections to be sure that integers are used
% when finding the middle of the shape)
if mod(length(int_x),2) == 1
    m = 0.5;
else
    m = 1;
end

% Get points for top and bottom rind sections
top_x = [ext_x(1:length(ext_x)/2+m),flip(int_x(1:length(ext_x)/2+m))];
top_y = [ext_y(1:length(ext_y)/2+m),flip(int_y(1:length(ext_y)/2+m))];
bot_x = [ext_x(length(ext_x)/2+m:length(ext_x)),flip(int_x(length(int_x)/2+m:length(int_x)))];
bot_y = [ext_y(length(ext_y)/2+m:length(ext_y)),flip(int_y(length(int_y)/2+m:length(int_y)))];

% calculate geometry and initial inertia values for top and bottom of cross section
[G_top,I_top,~] = polygeom(top_x,top_y);
[G_bot,I_bot,~] = polygeom(bot_x,bot_y);  

% Get inertia values about x axis
It = I_top(4);
Ib = I_bot(4);
It_temp = It;
Ib_temp = Ib;

% Get area and centroid vaues for top and bottom sections
At = G_top(1);
Ab = G_bot(1);
Ct = G_top(3);
Cb = G_bot(3);

It_i = It + At*(Ct)^2;
Ib_i = Ib + Ab*(Cb)^2;

% Determine which direction to shift center for first guess
if It_i - Ib_i > 0
    guess1 = 0;
    guess2 = 1;
elseif It_i - Ib_i < 0
    guess1 = -1;
    guess2 = 0;
end

% Initialize logic for iterations in 'while' loop
equal = 0;
n = 0;
guess_new = 0;

% Use a bisection method to narrow down 
while equal == 0 
    n = n + 1;
    % Check to see if inertias are within specified tolerance
    if abs(It_temp - Ib_temp) < tol
        equal = 1; % This is to end the while loop
    else
        % Use Parallel Axis Theorem to shift centroid. Use the same
        % halves with a shifted axis.
        % Determine how much to shift centroid
        guess_new = (guess2 + guess1)/2;
        
        % Calculate new inertia values based on new guess using the
        % parallel axis theorem
        It_temp = It + At*(Ct -guess_new)^2;
        Ib_temp = Ib + Ab*(-Cb + guess_new)^2;
            
        % Determine which side the inertia center is on
        if It_temp - Ib_temp > 0
            guess1 = guess_new;
        elseif It_temp - Ib_temp < 0
            guess2 = guess_new;
        end     
        
%         pause
        
    end           
end

% Record final centroid value
Inertia_Centroid_Y = guess_new;
Iuu = It_temp + Ib_temp;

end

