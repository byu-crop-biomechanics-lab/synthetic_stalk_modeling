load Section_slices_bottom_990.mat
load Ellipse_fits_bottom1.mat

for i = 1:length(indices)
    % Grab only the cross sections called out in indices vector
    section_ext = [ext_xDCSR(1,:,indices(i))', ext_yDCSR(1,:,indices(i))'];
    section_int = [int_xDCSR(1,:,indices(i))', int_yDCSR(1,:,indices(i))'];
    
    % Repeat the last points to close the loop
    section_ext = [section_ext; section_ext(1,:)];
    section_int = [section_int; section_int(1,:)];
    
    % Save base points as txt file
    S = size(section_ext);
    len = S(1);
    writespline(len,section_ext,i,'Real_ext')
    writespline(len,section_int,i,'Real_int')
    
    
end