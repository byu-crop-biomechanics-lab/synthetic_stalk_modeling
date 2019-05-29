clear;
close;
clear;

load Section_slices_bottom_990.mat
load Flip_sections.mat

ext_rhoDCSR = [];
int_rhoDCSR = [];
ext_xDCSR = ext_xDCSR;
ext_yDCSR = ext_yDCSR;
int_xDCSR = int_xDCSR;
int_yDCSR = int_yDCSR;
T = ext_tDCSR(:,1)';
N = size(ext_xDCSR,2);

flipped = []; % Hold the indices of the flipped sections
for i = 1:N
    if flip(i) ~= 0 
        flipped = [flipped; i];
        
        % Rotate and reorder external and internal points
        [~, ~, ~, ~, ~, ~, xp_ext, yp_ext, ~, ~] = reorder_V2(ext_xDCSR(:,i), ext_yDCSR(:,i), pi);
        [~, ~, ~, ~, ~, ~, xp_int, yp_int, ~, ~] = reorder_V2(int_xDCSR(:,i), int_yDCSR(:,i), pi);
        
        [~, ~, x_ext, y_ext, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_ext, yp_ext, 0);
        [~, ~, x_int, y_int, ~, ~, ~, ~, ~, ~] = reorder_V2(xp_int, yp_int, 0);
        
        % Redefine the appropriate row in the main XY data
        ext_xDCSR(:,i) = x_ext;
        ext_yDCSR(:,i) = y_ext;
        int_xDCSR(:,i) = x_int;
        int_yDCSR(:,i) = y_int;
    end
end

% % Plot the cross sections that should be flipped to verify
% for i = 1:length(flipped)
%     for j = 1:length(T)
%         plot(ext_xDCSR(1:j,flipped(i)),ext_yDCSR(1:j,flipped(i)));
%         hold on
%         plot(int_xDCSR(1:j,flipped(i)),int_yDCSR(1:j,flipped(i)));
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
        ext_rhoDCSR(i,j) = sqrt(ext_xDCSR(j,i)^2 + ext_yDCSR(j,i)^2);
        int_rhoDCSR(i,j) = sqrt(int_xDCSR(j,i)^2 + int_yDCSR(j,i)^2);
    end
end

% Transpose rho arrays so they are the same orientation as the other
% variables
ext_rhoDCSR = ext_rhoDCSR';
int_rhoDCSR = int_rhoDCSR';

% % Plot to check
% for i = 1:length(flipped)
%     polarplot(T,R_ext(flipped(i),:));
%     hold on
%     polarplot(T,R_int(flipped(i),:));
%     hold off
%     pause();
% end

% Save data as mat file
FolderName = pwd;
SaveName = 'Section_slices_bottom_990_FLIPPED.mat';
SaveFile = fullfile(FolderName, SaveName);
save(SaveFile,'ext_xDCSR','ext_yDCSR','int_xDCSR','int_yDCSR','ext_rhoDCSR','int_rhoDCSR','ext_tDCSR','int_tDCSR','avg_rind_thick');