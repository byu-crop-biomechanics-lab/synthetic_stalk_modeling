function remove_problem_ellipses(OriginalEllipseFits,problem_indices,GoodEllipseFits)
load(OriginalEllipseFits);

N = size(A,1);

% Remove the problem ellipses from the ellipse fit data
for i = 1:length(problem_indices)
    idx = problem_indices(i);
    if idx == 1
        A = [A(2:end,:)];
        AVG_RIND_T = [AVG_RIND_T(2:end)];
        B = [B(2:end)];
        DIFF_R_ext = [DIFF_R_ext(2:end,:)];
        DIFF_R_int = [DIFF_R_int(2:end,:)];
        ELLIPSE_CENTERS = [ELLIPSE_CENTERS(2:end,:)];
        ELLIPSE_R_ext = [ELLIPSE_R_ext(2:end,:)];
        ELLIPSE_R_int = [ELLIPSE_R_int(2:end,:)];
        ELLIPSE_T = [ELLIPSE_T(2:end,:)];
        ELLIPSE_XY = [ELLIPSE_XY(2:end,:,:)];
        R_ext = [R_ext(2:end,:)];
        R_int = [R_int(2:end,:)];
        
    elseif idx == N
        A = [A(1:N-1,:)];
        AVG_RIND_T = [AVG_RIND_T(1:N-1)];
        B = [B(1:N-1)];
        DIFF_R_ext = [DIFF_R_ext(1:N-1,:)];
        DIFF_R_int = [DIFF_R_int(1:N-1,:)];
        ELLIPSE_CENTERS = [ELLIPSE_CENTERS(1:N-1,:)];
        ELLIPSE_R_ext = [ELLIPSE_R_ext(1:N-1,:)];
        ELLIPSE_R_int = [ELLIPSE_R_int(1:N-1,:)];
        ELLIPSE_T = [ELLIPSE_T(1:N-1,:)];
        ELLIPSE_XY = [ELLIPSE_XY(1:N-1,:,:)];
        R_ext = [R_ext(1:N-1,:)];
        R_int = [R_int(1:N-1,:)];
        
    else
        A = [A(1:idx-1,:); A(idx+1:end,:)];
        AVG_RIND_T = [AVG_RIND_T(1:idx-1); AVG_RIND_T(idx+1:end)];
        B = [B(1:idx-1); B(idx+1:end)];
        DIFF_R_ext = [DIFF_R_ext(1:idx-1,:); DIFF_R_ext(idx+1:end,:)];
        DIFF_R_int = [DIFF_R_int(1:idx-1,:); DIFF_R_int(idx+1:end,:)];
        ELLIPSE_CENTERS = [ELLIPSE_CENTERS(1:idx-1,:); ELLIPSE_CENTERS(idx+1:end,:)];
        ELLIPSE_R_ext = [ELLIPSE_R_ext(1:idx-1,:); ELLIPSE_R_ext(idx+1:end,:)];
        ELLIPSE_R_int = [ELLIPSE_R_int(1:idx-1,:); ELLIPSE_R_int(idx+1:end,:)];
        ELLIPSE_T = [ELLIPSE_T(1:idx-1,:); ELLIPSE_T(idx+1:end,:)];
        ELLIPSE_XY = [ELLIPSE_XY(1:idx-1,:,:); ELLIPSE_XY(idx+1:end,:,:)];
        R_ext = [R_ext(1:idx-1,:); R_ext(idx+1:end,:)];
        R_int = [R_int(1:idx-1,:); R_int(idx+1:end,:)];
        
    end
    
    
end

% Save the final data in a new mat file
FolderName = pwd;
SaveFile       = fullfile(FolderName, GoodEllipseFits);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int',...
    'ELLIPSE_CENTERS','DIFF_R_ext','DIFF_R_int','R_ext','R_int','AVG_RIND_T','problem_indices');

end