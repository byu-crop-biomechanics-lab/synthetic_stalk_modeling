function [ERR_EXT,ERR_INT] = geometric_error(ChosenEllipseData,PCAData,NewGoodStalks,nNEPCs)
    
    hold off
    close all;
    set(0,'DefaultFigureWindowStyle','docked');

    load(PCAData,'ext_rhoPCAs','ext_rhocoeffs');
    load(ChosenEllipseData,'ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int','R_ext','R_int','AVG_RIND_T','B');
    load(NewGoodStalks,'newgoodstalknums');
    
    N = size(ELLIPSE_R_ext,1);
    npts = size(ELLIPSE_R_ext,2);
%     nNEPCs = 10;
    nrows = nNEPCs + 1;
    
    ERR_EXT = zeros(nrows,1);
    ERR_INT = zeros(nrows,1);
    
    
    % Work through NEPC cases:
    % Ellipse case
    [ERR_EXT(1),ERR_INT(1)] = get_total_scaled_error(N,npts,ELLIPSE_R_ext,ELLIPSE_R_int,R_ext,R_int,B);
    
    for row = 2:nrows
        NEPCnum = row - 1;

        NEPC_ext = zeros(size(ELLIPSE_R_ext));
        NEPC_int = zeros(size(ELLIPSE_R_ext));
        
        for i = 1:N
            [NEPC_ext(i,:),NEPC_int(i,:)] = section_from_PCA(newgoodstalknums,i,ELLIPSE_R_ext,ext_rhoPCAs,ext_rhocoeffs,AVG_RIND_T,NEPCnum);
        end
        
%         section = 6;
%         polarplot(ELLIPSE_T(section,:),NEPC_ext(section,:));
%         hold on
%         polarplot(ELLIPSE_T(section,:),R_ext(section,:));
%         hold off
%         pause();
%         close;
        
%         [ERR_EXT(row),ERR_INT(row)] = get_mean_error(N,npts,NEPC_ext,NEPC_int,R_ext,R_int);
        [ERR_EXT(row),ERR_INT(row)] = get_total_scaled_error(N,npts,NEPC_ext,NEPC_int,R_ext,R_int,B);
        
    end
    
%     % Normalize error
%     ERR_EXT_VAL = ERR_EXT(1);
%     ERR_INT_VAL = ERR_INT(1);
%     
%     for i = 1:nrows
%         ERR_EXT(i) = 100*ERR_EXT(i)/ERR_EXT_VAL;
%         ERR_INT(i) = 100*ERR_INT(i)/ERR_INT_VAL;
%     end


    % Plot results (normalized)
    cases = 0:1:nNEPCs;
%     plot(linspace(0,nNEPCs-1,nNEPCs+1),ERR_EXT);
    plot(cases,ERR_EXT);
    title('Average Geometric Error Progression');
    xlabel('Number of NEPCs Included');
    ylabel('Average Cumulative Radial Error');
end



%% Local functions
function [avg_err_ext,avg_err_int] = get_mean_error(N,npts,ELLIPSE_R_ext,ELLIPSE_R_int,R_ext,R_int)
    error_ext = zeros(N,npts);
    error_int = zeros(N,npts);
    
    for i = 1:N
        for j = 1:npts
            error_ext(i,j) = abs(ELLIPSE_R_ext(i,j) - R_ext(i,j));
            error_int(i,j) = abs(ELLIPSE_R_int(i,j) - R_int(i,j));            
        end
    end
    
%     plot(error_ext(i,:))
%     pause();


%     figure(1)
%     histogram(error_ext,8);
%     title('Exterior error');
%     
%     figure(2)
%     histogram(error_int,8);
%     title('Interior error');
%     
%     pause();
    
    % Take the average for each cross-section
    avg_err_ext = mean(sum(error_ext,2));
    avg_err_int = mean(sum(error_ext,2));

end


function [total_scaled_err_ext,total_scaled_err_int] = get_total_scaled_error(N,npts,ELLIPSE_R_ext,ELLIPSE_R_int,R_ext,R_int,B)
    error_ext = zeros(N,npts);
    error_int = zeros(N,npts);
    
    for i = 1:N
        for j = 1:npts
            error_ext(i,j) = abs(ELLIPSE_R_ext(i,j) - R_ext(i,j));
            error_int(i,j) = abs(ELLIPSE_R_int(i,j) - R_int(i,j));            
        end
    end
    
%     plot(error_ext(i,:))
%     pause();


%     figure(1)
%     histogram(error_ext,8);
%     title('Exterior error');
%     
%     figure(2)
%     histogram(error_int,8);
%     title('Interior error');
%     
%     pause();
    
    % Take the sum for each cross-section
    scaled_err_ext = sum(error_ext,2); % This is a vector
    scaled_err_int = sum(error_ext,2);
    
    for i = 1:N
        scaled_err_ext(i) = scaled_err_ext(i)/B(i);
        scaled_err_int(i) = scaled_err_int(i)/B(i);
    end
    
    total_scaled_err_ext = sum(scaled_err_ext);
    total_scaled_err_int = sum(scaled_err_int);

end


function [Rnew_ext,Rnew_int] = section_from_PCA(newgoodstalks,i,ELLIPSE_R_ext,ext_rhoPCAs,ext_rhocoeffs,AVG_RIND_T,numNEPCs)
    % THIS FUNCTION WORKS AS INTENDED WHEN USING THE PROPER DATA (LOAD PCA
    % DATA BEFORE ELLIPSE DATA SO ELLIPSE DATA OVERWRITES THE IMPORTANT
    % DATA IN THE CORRECT DIMENSIONS)
    
    NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
    
    for k = 1:numNEPCs
        % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
        NEPC_ext = NEPC_ext + ext_rhocoeffs(newgoodstalks(i),k)*ext_rhoPCAs(:,k)';
    end

    Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
    Rnew_int = Rnew_ext - AVG_RIND_T(i);

end