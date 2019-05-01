function ellipse_sensitivity(FileName)
% Take in a mat file with ellipse fit data (such as fivefits_top1.mat) and
% create new data points that represent a shift in either a, b, or t, as
% well as the original data points.

FolderName = pwd;
File       = fullfile(FolderName, FileName);
load(File,'A','B','ELLIPSE_XY','ELLIPSE_T','ELLIPSE_R','DIFF_ext_R','DIFF_int_R','AVG_RIND_T');

%% Save base XY points as CSV files
S = size(ELLIPSE_XY);
N = S(1);
for i = 1:N
    % Write exterior points
    ellipse_ext = [ELLIPSE_XY(i,:,1)', ELLIPSE_XY(i,:,2)'];
    
    % Write interior points (from polar data)
    ellipse_int = zeros(360,2);
    ELLIPSE_R_INT = ELLIPSE_R(i,:) - AVG_RIND_T(i);
    ELLIPSE_T_INT = ELLIPSE_T(i,:);
    
    for j = 1:360
        ellipse_int = convert_to_xy(ELLIPSE_R_INT, ELLIPSE_T_INT);
    end
    
    % Repeat the last points to close the loop
    ellipse_ext = [ellipse_ext; ellipse_ext(end,:)];
    ellipse_int = [ellipse_int; ellipse_int(end,:)];
    
    % Save base points as txt file
    S = size(ellipse_ext);
    len = S(1);
    writespline(len,ellipse_ext,i,'Base_ext')
    writespline(len,ellipse_int,i,'Base_int')
    
    
    
%     % Save base points as CSV
%     filename = sprintf('Base_ext%d.csv',i);
%     csvwrite(filename,ellipse_ext);
%     filename = sprintf('Base_int%d.csv',i);
%     csvwrite(filename,ellipse_int);
end

%%%%%%%%%%%%%%%%%% Sensitivity tests as CSV files %%%%%%%%%%%%%%%%%%%%%%%%%
%% CHANGING A
for i = 1:N
    % Adjust exterior points in polar
    Tnew = ELLIPSE_T(i,:);
    Aplus5_ext      = rpts(360,ELLIPSE_T(i,:),(1.05*A(i)),B(i));
    Aplus10_ext     = rpts(360,ELLIPSE_T(i,:),(1.10*A(i)),B(i));
    Aplus15_ext     = rpts(360,ELLIPSE_T(i,:),(1.15*A(i)),B(i));
    Aminus5_ext     = rpts(360,ELLIPSE_T(i,:),(0.95*A(i)),B(i));
    Aminus10_ext    = rpts(360,ELLIPSE_T(i,:),(0.90*A(i)),B(i));
    Aminus15_ext    = rpts(360,ELLIPSE_T(i,:),(0.85*A(i)),B(i));
    
    % Calculate the interior points
    Aplus5_int      = Aplus5_ext - AVG_RIND_T(i);
    Aplus10_int     = Aplus10_ext - AVG_RIND_T(i);
    Aplus15_int     = Aplus15_ext - AVG_RIND_T(i);
    Aminus5_int     = Aminus5_ext - AVG_RIND_T(i);
    Aminus10_int    = Aminus10_ext - AVG_RIND_T(i);
    Aminus15_int    = Aminus15_ext - AVG_RIND_T(i);
    
    % Convert to Cartesian
    Aplus5_xy_ext   = convert_to_xy(Aplus5_ext,Tnew);
    Aplus10_xy_ext  = convert_to_xy(Aplus10_ext,Tnew);
    Aplus15_xy_ext  = convert_to_xy(Aplus15_ext,Tnew);
    Aminus5_xy_ext  = convert_to_xy(Aminus5_ext,Tnew);
    Aminus10_xy_ext = convert_to_xy(Aminus10_ext,Tnew);
    Aminus15_xy_ext = convert_to_xy(Aminus15_ext,Tnew);
    
    Aplus5_xy_int   = convert_to_xy(Aplus5_int,Tnew);
    Aplus10_xy_int  = convert_to_xy(Aplus10_int,Tnew);
    Aplus15_xy_int  = convert_to_xy(Aplus15_int,Tnew);
    Aminus5_xy_int  = convert_to_xy(Aminus5_int,Tnew);
    Aminus10_xy_int = convert_to_xy(Aminus10_int,Tnew);
    Aminus15_xy_int = convert_to_xy(Aminus15_int,Tnew);
    
    % Repeat the last points to close the loop
    Aplus5_xy_ext   = [Aplus5_xy_ext; Aplus5_xy_ext(end,:)];
    Aplus10_xy_ext  = [Aplus10_xy_ext; Aplus10_xy_ext(end,:)];
    Aplus15_xy_ext  = [Aplus15_xy_ext; Aplus15_xy_ext(end,:)];
    Aminus5_xy_ext  = [Aminus5_xy_ext; Aminus5_xy_ext(end,:)];
    Aminus10_xy_ext = [Aminus10_xy_ext; Aminus10_xy_ext(end,:)];
    Aminus15_xy_ext = [Aminus15_xy_ext; Aminus15_xy_ext(end,:)];
    
    Aplus5_xy_int   = [Aplus5_xy_int; Aplus5_xy_int(end,:)];
    Aplus10_xy_int  = [Aplus10_xy_int; Aplus10_xy_int(end,:)];
    Aplus15_xy_int  = [Aplus15_xy_int; Aplus15_xy_int(end,:)];
    Aminus5_xy_int  = [Aminus5_xy_int; Aminus5_xy_int(end,:)];
    Aminus10_xy_int = [Aminus10_xy_int; Aminus10_xy_int(end,:)];
    Aminus15_xy_int = [Aminus15_xy_int; Aminus15_xy_int(end,:)];
    
    % Save base points as txt file
    S = size(Aplus5_xy_ext);
    len = S(1);
    writespline(len,Aplus5_xy_ext,i,'Aplus5_ext');
    writespline(len,Aplus10_xy_ext,i,'Aplus10_ext');
    writespline(len,Aplus15_xy_ext,i,'Aplus15_ext');
    writespline(len,Aminus5_xy_ext,i,'Aminus5_ext');
    writespline(len,Aminus10_xy_ext,i,'Aminus10_ext');
    writespline(len,Aminus15_xy_ext,i,'Aminus15_ext');
    
    writespline(len,Aplus5_xy_int,i,'Aplus5_int');
    writespline(len,Aplus10_xy_int,i,'Aplus10_int');
    writespline(len,Aplus15_xy_int,i,'Aplus15_int');
    writespline(len,Aminus5_xy_int,i,'Aminus5_int');
    writespline(len,Aminus10_xy_int,i,'Aminus10_int');
    writespline(len,Aminus15_xy_int,i,'Aminus15_int');
    
%     % Save new xy points as CSV files
%     filename = sprintf('Aplus5_ext%d.csv',i);
%     csvwrite(filename,Aplus5_xy_ext);
%     filename = sprintf('Aplus10_ext%d.csv',i);
%     csvwrite(filename,Aplus10_xy_ext);
%     filename = sprintf('Aplus15_ext%d.csv',i);
%     csvwrite(filename,Aplus15_xy_ext);
%     filename = sprintf('Aminus5_ext%d.csv',i);
%     csvwrite(filename,Aminus5_xy_ext);
%     filename = sprintf('Aminus10_ext%d.csv',i);
%     csvwrite(filename,Aminus10_xy_ext);
%     filename = sprintf('Aminus15_ext%d.csv',i);
%     csvwrite(filename,Aminus15_xy_ext);
%     
%     filename = sprintf('Aplus5_int%d.csv',i);
%     csvwrite(filename,Aplus5_xy_int);
%     filename = sprintf('Aplus10_int%d.csv',i);
%     csvwrite(filename,Aplus10_xy_int);
%     filename = sprintf('Aplus15_int%d.csv',i);
%     csvwrite(filename,Aplus15_xy_int);
%     filename = sprintf('Aminus5_int%d.csv',i);
%     csvwrite(filename,Aminus5_xy_int);
%     filename = sprintf('Aminus10_int%d.csv',i);
%     csvwrite(filename,Aminus10_xy_int);
%     filename = sprintf('Aminus15_int%d.csv',i);
%     csvwrite(filename,Aminus15_xy_int);
end

%% CHANGING B
for i = 1:N
    % Adjust exterior points in polar
    Tnew = ELLIPSE_T(i,:);
    Bplus5_ext      = rpts(360,ELLIPSE_T(i,:),(A(i)),1.05*B(i));
    Bplus10_ext     = rpts(360,ELLIPSE_T(i,:),(A(i)),1.10*B(i));
    Bplus15_ext     = rpts(360,ELLIPSE_T(i,:),(A(i)),1.15*B(i));
    Bminus5_ext     = rpts(360,ELLIPSE_T(i,:),(A(i)),0.95*B(i));
    Bminus10_ext    = rpts(360,ELLIPSE_T(i,:),(A(i)),0.90*B(i));
    Bminus15_ext    = rpts(360,ELLIPSE_T(i,:),(A(i)),0.85*B(i));
    
    % Calculate the interior points
    Bplus5_int      = Bplus5_ext - AVG_RIND_T(i);
    Bplus10_int     = Bplus10_ext - AVG_RIND_T(i);
    Bplus15_int     = Bplus15_ext - AVG_RIND_T(i);
    Bminus5_int     = Bminus5_ext - AVG_RIND_T(i);
    Bminus10_int    = Bminus10_ext - AVG_RIND_T(i);
    Bminus15_int    = Bminus15_ext - AVG_RIND_T(i);
    
    % Convert to Cartesian
    Bplus5_xy_ext   = convert_to_xy(Bplus5_ext,Tnew);
    Bplus10_xy_ext  = convert_to_xy(Bplus10_ext,Tnew);
    Bplus15_xy_ext  = convert_to_xy(Bplus15_ext,Tnew);
    Bminus5_xy_ext  = convert_to_xy(Bminus5_ext,Tnew);
    Bminus10_xy_ext = convert_to_xy(Bminus10_ext,Tnew);
    Bminus15_xy_ext = convert_to_xy(Bminus15_ext,Tnew);
    
    Bplus5_xy_int   = convert_to_xy(Bplus5_int,Tnew);
    Bplus10_xy_int  = convert_to_xy(Bplus10_int,Tnew);
    Bplus15_xy_int  = convert_to_xy(Bplus15_int,Tnew);
    Bminus5_xy_int  = convert_to_xy(Bminus5_int,Tnew);
    Bminus10_xy_int = convert_to_xy(Bminus10_int,Tnew);
    Bminus15_xy_int = convert_to_xy(Bminus15_int,Tnew);
    
    % Repeat the last points to close the loop
    Bplus5_xy_ext   = [Bplus5_xy_ext; Bplus5_xy_ext(end,:)];
    Bplus10_xy_ext  = [Bplus10_xy_ext; Bplus10_xy_ext(end,:)];
    Bplus15_xy_ext  = [Bplus15_xy_ext; Bplus15_xy_ext(end,:)];
    Bminus5_xy_ext  = [Bminus5_xy_ext; Bminus5_xy_ext(end,:)];
    Bminus10_xy_ext = [Bminus10_xy_ext; Bminus10_xy_ext(end,:)];
    Bminus15_xy_ext = [Bminus15_xy_ext; Bminus15_xy_ext(end,:)];
    
    Bplus5_xy_int   = [Bplus5_xy_int; Bplus5_xy_int(end,:)];
    Bplus10_xy_int  = [Bplus10_xy_int; Bplus10_xy_int(end,:)];
    Bplus15_xy_int  = [Bplus15_xy_int; Bplus15_xy_int(end,:)];
    Bminus5_xy_int  = [Bminus5_xy_int; Bminus5_xy_int(end,:)];
    Bminus10_xy_int = [Bminus10_xy_int; Bminus10_xy_int(end,:)];
    Bminus15_xy_int = [Bminus15_xy_int; Bminus15_xy_int(end,:)];
    
    % Save base points as txt file
    S = size(Bplus5_xy_ext);
    len = S(1);
    writespline(len,Bplus5_xy_ext,i,'Bplus5_ext');
    writespline(len,Bplus10_xy_ext,i,'Bplus10_ext');
    writespline(len,Bplus15_xy_ext,i,'Bplus15_ext');
    writespline(len,Bminus5_xy_ext,i,'Bminus5_ext');
    writespline(len,Bminus10_xy_ext,i,'Bminus10_ext');
    writespline(len,Bminus15_xy_ext,i,'Bminus15_ext');
    
    writespline(len,Bplus5_xy_int,i,'Bplus5_int');
    writespline(len,Bplus10_xy_int,i,'Bplus10_int');
    writespline(len,Bplus15_xy_int,i,'Bplus15_int');
    writespline(len,Bminus5_xy_int,i,'Bminus5_int');
    writespline(len,Bminus10_xy_int,i,'Bminus10_int');
    writespline(len,Bminus15_xy_int,i,'Bminus15_int');
    
%     % Save new xy points as CSV files
%     filename = sprintf('Bplus5_ext%d.csv',i);
%     csvwrite(filename,Bplus5_xy_ext);
%     filename = sprintf('Bplus10_ext%d.csv',i);
%     csvwrite(filename,Bplus10_xy_ext);
%     filename = sprintf('Bplus15_ext%d.csv',i);
%     csvwrite(filename,Bplus15_xy_ext);
%     filename = sprintf('Bminus5_ext%d.csv',i);
%     csvwrite(filename,Bminus5_xy_ext);
%     filename = sprintf('Bminus10_ext%d.csv',i);
%     csvwrite(filename,Bminus10_xy_ext);
%     filename = sprintf('Bminus15_ext%d.csv',i);
%     csvwrite(filename,Bminus15_xy_ext);
%     
%     filename = sprintf('Bplus5_int%d.csv',i);
%     csvwrite(filename,Bplus5_xy_int);
%     filename = sprintf('Bplus10_int%d.csv',i);
%     csvwrite(filename,Bplus10_xy_int);
%     filename = sprintf('Bplus15_int%d.csv',i);
%     csvwrite(filename,Bplus15_xy_int);
%     filename = sprintf('Bminus5_int%d.csv',i);
%     csvwrite(filename,Bminus5_xy_int);
%     filename = sprintf('Bminus10_int%d.csv',i);
%     csvwrite(filename,Bminus10_xy_int);
%     filename = sprintf('Bminus15_int%d.csv',i);
%     csvwrite(filename,Bminus15_xy_int);
end

%% CHANGING T
for i = 1:N
    Tnew = ELLIPSE_T(i,:);
    base_ext = rpts(360,ELLIPSE_T(i,:),(A(i)),B(i));
    
    % Calculate the interior points
    Tplus5_int      = base_ext - 1.05*AVG_RIND_T(i);
    Tplus10_int     = base_ext - 1.10*AVG_RIND_T(i);
    Tplus15_int     = base_ext - 1.15*AVG_RIND_T(i);
    Tminus5_int     = base_ext - 0.95*AVG_RIND_T(i);
    Tminus10_int    = base_ext - 0.90*AVG_RIND_T(i);
    Tminus15_int    = base_ext - 0.85*AVG_RIND_T(i);
    
    % Convert to Cartesian
    base_xy_ext     = convert_to_xy(base_ext,Tnew);
    Tplus5_xy_int   = convert_to_xy(Tplus5_int,Tnew);
    Tplus10_xy_int  = convert_to_xy(Tplus10_int,Tnew);
    Tplus15_xy_int  = convert_to_xy(Tplus15_int,Tnew);
    Tminus5_xy_int  = convert_to_xy(Tminus5_int,Tnew);
    Tminus10_xy_int = convert_to_xy(Tminus10_int,Tnew);
    Tminus15_xy_int = convert_to_xy(Tminus15_int,Tnew);
    
    % Repeat the last points to close the loop
    base_xy_ext   = [base_xy_ext; base_xy_ext(end,:)];
    Tplus5_xy_int   = [Tplus5_xy_int; Tplus5_xy_int(end,:)];
    Tplus10_xy_int  = [Tplus10_xy_int; Tplus10_xy_int(end,:)];
    Tplus15_xy_int  = [Tplus15_xy_int; Tplus15_xy_int(end,:)];
    Tminus5_xy_int  = [Tminus5_xy_int; Tminus5_xy_int(end,:)];
    Tminus10_xy_int = [Tminus10_xy_int; Tminus10_xy_int(end,:)];
    Tminus15_xy_int = [Tminus15_xy_int; Tminus15_xy_int(end,:)];
    
    % Save base points as txt file
    S = size(base_xy_ext);
    len = S(1);
    writespline(len,base_xy_ext,i,'base_ext');
    writespline(len,Tplus5_xy_int,i,'Tplus5_int');
    writespline(len,Tplus10_xy_int,i,'Tplus10_int');
    writespline(len,Tplus15_xy_int,i,'Tplus15_int');
    writespline(len,Tminus5_xy_int,i,'Tminus5_int');
    writespline(len,Tminus10_xy_int,i,'Tminus10_int');
    writespline(len,Tminus15_xy_int,i,'Tminus15_int');
    
%     % Save new xy points as CSV files
%     filename = sprintf('base_xy_ext%d.csv',i);
%     csvwrite(filename,base_xy_ext);
%     filename = sprintf('Tplus5_xy_int%d.csv',i);
%     csvwrite(filename,Tplus5_xy_int);
%     filename = sprintf('Tplus10_xy_int%d.csv',i);
%     csvwrite(filename,Tplus10_xy_int);
%     filename = sprintf('Tplus15_xy_int%d.csv',i);
%     csvwrite(filename,Tplus15_xy_int);
%     filename = sprintf('Tminus5_xy_int%d.csv',i);
%     csvwrite(filename,Tminus5_xy_int);
%     filename = sprintf('Tminus10_xy_int%d.csv',i);
%     csvwrite(filename,Tminus10_xy_int);
%     filename = sprintf('Tminus15_xy_int%d.csv',i);
%     csvwrite(filename,Tminus15_xy_int);

end

end



function [r] = rpts(N,theta,dmaj,dmin)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2);
    end
end

function [xy_columns] = convert_to_xy(R,theta)
    N = length(theta);
    xy_columns = zeros(N,2);
    for i = 1:N
        xy_columns(i,1) = R(i)*cos(theta(i));
        xy_columns(i,2) = R(i)*sin(theta(i));
    end
end

function writespline(len,data,k,outputName)
    %define empty spline and number of x-y points
    spline = '';

    %run through 1-column arrays of the x and y data points for the spline, and add to the end of the string with the correct formatting
    for i = 1:len 
        spline = strcat(spline,'(',num2str(data(i,1)),', ',num2str(data(i,2)),'), '); 
    end

    %write string to file using the specimen name
    I = sprintf('%d',k);
    fid = fopen(strcat(outputName,I,'.txt'),'wt');
    fprintf(fid, spline);
    fclose(fid);
end