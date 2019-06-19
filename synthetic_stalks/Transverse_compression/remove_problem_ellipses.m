function remove_problem_ellipses(OriginalEllipseFits,problem_indices,GoodEllipseFits)
load(OriginalEllipseFits);

N = size(A,1);

% Remove the problem ellipses from the ellipse fit data
A(problem_indices) = [];
AVG_RIND_T(problem_indices) = [];
B(problem_indices) = [];
DIFF_R_ext(problem_indices,:) = [];
DIFF_R_int(problem_indices,:) = [];
ELLIPSE_CENTERS(problem_indices,:) = [];
ELLIPSE_R_ext(problem_indices,:) = [];
ELLIPSE_R_int(problem_indices,:) = [];
ELLIPSE_T(problem_indices,:) = [];
ELLIPSE_XY(problem_indices,:,:) = [];
R_ext(problem_indices,:) = [];
R_int(problem_indices,:) = [];

% Save the final data in a new mat file
FolderName = pwd;
SaveFile       = fullfile(FolderName, GoodEllipseFits);
save(SaveFile,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int',...
    'ELLIPSE_CENTERS','DIFF_R_ext','DIFF_R_int','R_ext','R_int','AVG_RIND_T','problem_indices');

end