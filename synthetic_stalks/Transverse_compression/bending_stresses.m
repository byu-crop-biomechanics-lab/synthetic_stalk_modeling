% bending_stresses.m: Calculate the bending stresses in a stalk
% cross-section using its area and material properties. Uses polygeom.m
% from Matlab File Exchange
% 
% Author: Ryan Larson
% Date: 6/11/2019
close; clear;

load NEPCs_bottom_987.mat
load Ellipse_fits_bottom_987.mat

% Material properties
Er = 14672.916666667;
Ep = 1467.2916666667;

% Moment load
M = 100000; % N*mm

% Create arrays for holding rind and pith results
rind_stresses = NaN(50,11);
rind_Ixx = NaN(50,11);
rind_A = NaN(50,11);
rind_centroid = NaN(50,11,2);

pith_stresses = NaN(50,11);
pith_Ixx = NaN(50,11);
pith_A = NaN(50,11);
pith_centroid = NaN(50,11,2);

% Iterate through the cross-sections
for i = 1:50

    % Iterate through the cases
    % Real (Case 0)
    
% Convert data to Cartesian coordinates (read in as row vectors)
if size(R_ext,1) > 1
    X_ext = R_ext(i,:).*cos(ELLIPSE_T(i,:));
    Y_ext = R_ext(i,:).*sin(ELLIPSE_T(i,:));
    X_int = R_int(i,:).*cos(ELLIPSE_T(i,:));
    Y_int = R_int(i,:).*sin(ELLIPSE_T(i,:));
% else
%     X_ext = R_ext(1,:).*cos(T(1,:));
%     Y_ext = R_ext(1,:).*sin(T(1,:));
%     X_int = R_int(1,:).*cos(T(1,:));
%     Y_int = R_int(1,:).*sin(T(1,:));
end

X_ext = [X_ext,X_ext(1)];
Y_ext = [Y_ext,Y_ext(1)];
X_int = [X_int,X_int(1)];
Y_int = [Y_int,Y_int(1)];

X_int_CW = flip(X_int);
Y_int_CW = flip(Y_int);

% Pith geometry data
[ int_geom, int_iner, int_cpmo ] = polygeom(X_int,Y_int);
Ixx_p = int_iner(1)
Ap = int_geom(1)
Xc_p = int_geom(2)
Yc_p = int_geom(3)

% Rind geometry data
X_ext_temp = [X_ext, X_int_CW, X_ext(1)];
Y_ext_temp = [Y_ext, Y_int_CW, Y_ext(1)];
[ ext_geom, ext_iner, ext_cpmo ] = polygeom(X_ext_temp,Y_ext_temp);
Ixx_r = ext_iner(1)
Ar = ext_geom(1)
Xc_r = ext_geom(2)
Yc_r = ext_geom(3)

%% Determine the neutral axis based on geometry
% Distance in y from minimum y point to centers of mass for rind and pith
min_y = min(Y_ext)
dr = abs(Yc_r - min_y)
dp = abs(Yc_p - min_y)

% Calculate h (distance from lowest y point to the neutral axis) from
% areas, material properties, and d distances
h = (Er*dr*Ar + Ep*dp*Ap)/(Er*Ar + Ep*Ap)

% Max stress in the pith
[maxp, index] = max(abs(Y_int));
maxp = maxp*sign(Y_int(index));
sigma_p = (-M*maxp*Ep)/(Ep*Ixx_p + Er*Ixx_r)

% Max stress in the rind
[maxr, index] = max(abs(Y_ext));
maxr = maxr*sign(Y_ext(index));
sigma_r = (-M*maxr*Er)/(Ep*Ixx_p + Er*Ixx_r)
