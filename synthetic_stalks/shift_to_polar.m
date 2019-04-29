function shift_to_polar(FileName,angle,SaveName)
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
    
    R = [cosd(angle), -sind(angle);
         sind(angle), cosd(angle)];
    
    % Rotate exterior points
    for i = 1:length(ext_X)
        x = ext_X(i);
        y = ext_Y(i);
        
        PT = [x;y];
        
        PTrot = R*PT;
        
        ext_X(i) = PTrot(1);
        ext_Y(i) = PTrot(2);
    end
    
    % Rotate interior points
    for i = 1:length(int_X)
        x = int_X(i);
        y = int_Y(i);
        
        PT = [x;y];
        
        PTrot = R*PT;
        
        int_X(i) = PTrot(1);
        int_Y(i) = PTrot(2);
    end
    
    alpha = 0;
    [~, ~, ~, ~, ~, ~, ~, ~, ext_R, ext_T] = reorder_V2(ext_X, ext_Y, alpha);
    [~, ~, ~, ~, ~, ~, ~, ~, int_R, int_T] = reorder_V2(int_X, int_Y, alpha);
    
    % Check the points for any repeats in theta. If there is a repeat, then
    % shift the second instance by 0.001.
%     disp('EXTERIOR');
    for i = 2:length(ext_T)
        if ext_T(i) == ext_T(i-1)
            ext_T(i) = ext_T(i) + 0.0001;
%             disp('repeated value')
        end    
    end
%     disp('INTERIOR');

    for i = 2:length(int_T)
        if int_T(i) == int_T(i-1)
            int_T(i) = int_T(i) + 0.0001;
%             disp('repeated value')
        end    
    end

    % Downsample the exterior and interior points
    [ext_R_sample, ext_T_sample, ~, ~] = downsampler(ext_R, ext_T, 0, 0, 360);
    [int_R_sample, int_T_sample, ~, ~] = downsampler(int_R, int_T, 0, 0, 360);

    % Close the loop by repeating the first point in all vectors
    ext_T_sample = [ext_T_sample; ext_T_sample(1)];
    ext_R_sample = [ext_R_sample; ext_R_sample(1)];
    int_T_sample = [int_T_sample; int_T_sample(1)];
    int_R_sample = [int_R_sample; int_R_sample(1)];

    ext_T = ext_T_sample;
    ext_R = ext_R_sample;
    int_T = int_T_sample;
    int_R = int_R_sample;

    polarplot(ext_T,ext_R);
    hold on
    polarplot(int_T,int_R);

    SaveFile       = fullfile(FolderName, SaveName);
    save(SaveFile,'ext_T','ext_R','int_T','int_R');
end