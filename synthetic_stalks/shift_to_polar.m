function [ext_T,ext_R,int_T,int_R] = shift_to_polar(FileName,SaveName)
% Take in boundary data from the .mat file specified by filename, shift it
% to be centered at the origin, and rotate it by the angle specified in
% Chris's spreadsheet. Then convert to polar coordinates and save as a new
% .mat file.
    
    FolderName = pwd;
    File       = fullfile(FolderName, FileName);
    load(File,'ext_X','ext_Y','int_X','int_Y','xbar','ybar');   % not: load('File')

    ext_X = ext_X - xbar;
    ext_Y = ext_Y - ybar;
    int_X = int_X - xbar;
    int_Y = int_Y - ybar;
    
%     R = [cosd(angle), -sind(angle);
%          sind(angle), cosd(angle)];
%     
%     % Rotate exterior points
%     for i = 1:length(ext_X)
%         x = ext_X(i);
%         y = ext_Y(i);
%         
%         PT = [x;y];
%         
%         PTrot = R*PT;
%         
%         ext_X(i) = PTrot(1);
%         ext_Y(i) = PTrot(2);
%     end
%     
%     % Rotate interior points
%     for i = 1:length(int_X)
%         x = int_X(i);
%         y = int_Y(i);
%         
%         PT = [x;y];
%         
%         PTrot = R*PT;
%         
%         int_X(i) = PTrot(1);
%         int_Y(i) = PTrot(2);
%     end
    
    [ext_T,ext_R] = cart2pol(ext_X,ext_Y);
    [int_T,int_R] = cart2pol(int_X,int_Y);
    
    polarplot(ext_T,ext_R);
    hold on
    polarplot(int_T,int_R);
    
    SaveFile       = fullfile(FolderName, SaveName);
    save(SaveFile,'ext_T','ext_R','int_T','int_R');
end