 function [] = generate_curves(stalk, thresh, zsample, filename_ext, filename_int)

%=======================================================================================================
% INPUTS 
%   stalk -         a 3D grayscale data file (i.e., "stalk" variable from a single stalk .mat file)
%   threshold -     greyscale threshold value used for exterior boundary
%                   deterection (currently using 15000)
%   zsample -       how many slices to skip between samples (e.g. "2" means
%                   take every other slice)
%   filename_ext -  the name of the file, including the .ply extension (e.g.
%                   'example.ply'), for the external curve
%   filename_in -   the name of the file, including the .ply extension (e.g.
%                   'example.ply'), for the internal curve
%   example -       generate_curves(stalk, 15000,10,'ext.ply','int.ply')


struct_curve_ext = [];
struct_curve_int = [];
Nslices = size(stalk,3);                    % number of cross-sectional slices in this data structure
plotting = 0;                               % no plotting option


        
for i = 100:zsample:Nslices-100                   % Loop through each cross-sectional slice
    I = stalk(:,:,i);                       % extract one cross-sectional image, "I"
    
    [ ext_X, ext_Y, ~, ~] = exterior_boundaries_V4( I, thresh, plotting );

    local_curve_ext = [];
    local_curve_ext(:,1) = ext_X;
    local_curve_ext(:,2) = ext_Y;
    local_curve_ext(:,3) = i;
    struct_curve_ext = [struct_curve_ext;local_curve_ext];
    
    [~, int_X, int_Y] = avg_rind_thickness_normal_method(I, ext_X, ext_Y, plotting);
    
    local_curve_int = [];
    local_curve_int(:,1) = int_X;
    local_curve_int(:,2) = int_Y;
    local_curve_int(:,3) = i;
    struct_curve_int = [struct_curve_int;local_curve_int];
    
end

% write the files

ptCloud_ext = pointCloud(struct_curve_ext);
pcwrite(ptCloud_ext,filename_ext);
ptCloud_int = pointCloud(struct_curve_int);
pcwrite(ptCloud_int,filename_int);
