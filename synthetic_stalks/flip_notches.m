clear;
close;

load Section_slices_bottom_990.mat
load Flip_sections.mat

R_ext = squeeze(ext_rhoDCSR);
R_int = squeeze(int_rhoDCSR);
X_ext = squeeze(ext_xDCSR);
Y_ext = squeeze(ext_yDCSR);
X_int = squeeze(int_xDCSR);
Y_int = squeeze(int_yDCSR);
T = tDCSR';
N = length(T);

flipped = []; % Hold the indices of the flipped sections
for i = 1:N
    if flip(i) ~= 0
        flipped = [flipped; i];
        
        % Rotate and reorder external and internal points
        [~, ~, ~, ~, ~, ~, xp_ext, yp_ext, ~, ~] = reorder_V2(X_ext(:,i), Y_ext(:,i), pi);
        [~, ~, ~, ~, ~, ~, xp_int, yp_int, ~, ~] = reorder_V2(X_int(:,i), Y_int(:,i), pi);
        
        [~, ~, x_ext, y_ext, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_ext, yp_ext, 0);
        [~, ~, x_int, y_int, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_int, yp_int, 0);
        
        % Redefine the appropriate row in the main XY data
        X_ext(:,i) = x_ext;
        Y_ext(:,i) = y_ext;
        X_int(:,i) = x_int;
        Y_int(:,i) = y_int;
    end
end

% % Plot the cross sections that should be flipped to verify
% for i = 1:length(flipped)
%     for j = 1:length(T)
%         plot(X_ext(1:j,flipped(i)),Y_ext(1:j,flipped(i)));
%         hold on
%         plot(X_int(1:j,flipped(i)),Y_int(1:j,flipped(i)));
%         axis equal
%         hold off
%         pause(0.01);
%         
%     end
%     pause();
% end


% Convert data from Cartesian to polar
for i = 1:N
    for j = 1:length(T)
        R_ext(i,j) = sqrt(X_ext(i,j)^2 + Y_ext(i,j)^2);
        R_int(i,j) = sqrt(X_int(i,j)^2 + Y_int(i,j)^2);
    end
end

% % Plot to check
% for i = 1:length(flipped)
%     polarplot(T,R_ext(:,flipped(i)));
%     hold on
%     polarplot(T,R_int(:,flipped(i)));
%     hold off
%     pause();
% end

% Save data as mat file
FolderName = pwd;
SaveName = 'Section_slices_bottom_990_FLIPPED.mat';
SaveFile = fullfile(FolderName, SaveName);
save(SaveFile,'X_ext','Y_ext','X_int','Y_int','R_ext','R_int','T','avg_rind_thickness');