function [Ix_rind,Ix_pith] = moment_inertia(theta,int_r,ext_r)
% FILENAME: moment_inertia.m
% AUTHOR: Michael Ottesen
% DATE: 5/6/2020
%
% PURPOSE: Calcualte area moment of inertia for any stalk cross section.
% 
% INPUTS:
%       theta: Theta values for interior and exterior points. These are in
%       radians, taken from results from AllSlicesPCA
%
%       int_r: Internal radius points that define the cross section that 
%       correlate with theta
%
%       ext_r: External radius points that define the cross section that
%       correlate with theta.
%       
% OUTPUTS:
%       Ix_rind: real value representing the area moment of inertia of the 
%       rind about the x-x axis through the centroid of the cross section.
%
%       Ix_pith: real value representing the area moment of inertia of the 
%       pith about the x-x axis through the centroid of the cross section.
%
% NOTES: 
%       - polygeom (check on file exchange)
% 
% PSEUDO-CODE:
%   
% 
% VERSION HISTORY:
% V1 - Used different rotation method, but gave wrong answers
%     Ix_rec = (1/12)*b_rec*h_rec^3 + A(i)*(Cx(i)^2);
%     Iy_rec = (1/12)*b_rec^3*h_rec + A(i)*(Cy(i)^2);
%     Ixy_rec = A(i)*Cx(i)*Cy(i);
%     Ix_rotated = (Ix_rec+Iy_rec)/2 + ((Ix_rec-Iy_rec)/2)*cos(-2*theta(i)) - Ixy_rec*sin(-2*theta(i));
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------

% Determine number of segments
N = length(theta);

% Convert polar coordinates to cartesian for analysis
[int_x,int_y] = pol2cart(theta,int_r);
[ext_x,ext_y] = pol2cart(theta,ext_r);

% Add point at the end to be a repeat point of the first one
int_x = [int_x,int_x(1)];
int_y = [int_y,int_y(1)];
ext_x = [ext_x,ext_x(1)];
ext_y = [ext_y,ext_y(1)]; 

% Initialize Ix values
Ix_rind = 0; % Initialize Ix_rind
Ix_pith = 0; % Initialize Ix_pith

for i = 1:N
    %========================= Inertia of rind ==========================
    % Get points for quadrilateral
    Poly_x = [ext_x(i), ext_x(i+1), int_x(i+1), int_x(i)];
    Poly_y = [ext_y(i), ext_y(i+1), int_y(i+1), int_y(i)];
    
    % Find area and centroid of quadrilateral
    A(i) = polyarea(Poly_x,Poly_y);
    polyin = polyshape(Poly_x,Poly_y);
    [Cx(i), Cy(i)] = centroid(polyin);
    
    % Get base and height values of small rectangle using pythag theorem
    b1 = sqrt((ext_x(i)-int_x(i))^2 + (ext_y(i)-int_y(i))^2);
    b2 = sqrt((ext_x(i+1)-int_x(i+1))^2 + (ext_y(i+1)-int_y(i+1))^2);
    b_rec = (b1 + b2)/2;
    
    h1 = sqrt((ext_x(i)-ext_x(i+1))^2 + (ext_y(i)-ext_y(i+1))^2);
    h2 = sqrt((int_x(i)-int_x(i+1))^2 + (int_y(i)-int_y(i+1))^2);
    h_rec = (h1 + h2)/2;
    
    % Calculate moment of inertia of approximate rectangle. Then parallel
    % axis theorem is used to combine into total rind moment of inertia
    Ix_rotated = (1/12)*b_rec*h_rec*((h_rec^2)*cos(theta(i))^2 + (b_rec^2)*sin(theta(i))^2);
    Ix_PrinAxis = A(i)*(Cy(i)^2);
    Ix_segment = Ix_rotated + Ix_PrinAxis;
    Ix_rind = Ix_rind + Ix_segment;

    
    %========================= Inertia of pith ==========================
    % Cut pith into skinny traingles and find inertia values for each about
    % the x axis (they are rotated to be about same axis)
    
    % Get points for triangle
    Poly_x2 = [0, int_x(i), int_x(i+1)];
    Poly_y2 = [0, int_y(i), int_y(i+1)];
    
    % Find area and centroid of triangle
    A2(i) = polyarea(Poly_x2,Poly_y2);
    polyin2 = polyshape(Poly_x2,Poly_y2);
    [Cx2(i), Cy2(i)] = centroid(polyin2);
    
    % Find eaquation of the line from origin to first point
    coefficients = polyfit([0, int_x(i)], [0, int_y(i)], 1);
    a = -coefficients(1);
    b = 1;
    c = -coefficients(2);
    
    % find length of line, and distance from next point to line
    b_tri = sqrt(int_x(i)^2 + int_y(i)^2);
    h_tri = abs(a*int_x(i+1) + b*int_y(i+1) + c)/sqrt(a^2 + b^2); % point to line
    
    % Calculate moment of inertia of rotated triangle about global x axis
    Ix_rotated2 = (1/36)*b_tri*h_tri*((h_tri^2)*cos(theta(i))^2 + (b_tri^2)*sin(theta(i))^2);
    Ix_PrinAxis2 = A2(i)*(Cy2(i)^2);
    Ix_segment2 = Ix_rotated2 + Ix_PrinAxis2;
    Ix_pith = Ix_pith + Ix_segment2;

end
end
