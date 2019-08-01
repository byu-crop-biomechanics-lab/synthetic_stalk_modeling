function [ERR_EXT,ERR_INT] = geometric_error(ChosenEllipseData,PCAData)
    
    hold off
    close all;
    set(0,'DefaultFigureWindowStyle','docked');

    load(ChosenEllipseData,'ELLIPSE_T','ELLIPSE_R_ext','ELLIPSE_R_int','R_ext','R_int','AVG_RIND_T');
    load(PCAData,'ext_rhoPCAs','ext_rhocoeffs');
    
    N = size(ELLIPSE_R_ext,1);
    npts = size(ELLIPSE_R_ext,2);
    
    ERR_EXT = zeros(6,1);
    ERR_INT = zeros(6,1);
    
    
    
    % Work through NEPC cases:
    % Ellipse case
    [ERR_EXT(1),ERR_INT(1)] = get_mean_error(N,npts,ELLIPSE_R_ext,ELLIPSE_R_int,R_ext,R_int);
    
    for row = 2:6
        numNEPCs = row - 1;

        for i = 1:N
            NEPC_ext = zeros(size(ELLIPSE_R_ext));
            NEPC_int = zeros(size(ELLIPSE_R_ext));
            [NEPC_ext(i,:),NEPC_int(i,:)] = section_from_PCA(i,ELLIPSE_R_ext,ext_rhoPCAs,ext_rhocoeffs,AVG_RIND_T,numNEPCs);
        end
        
%         section = 1;
%         polarplot(ELLIPSE_T(section,:),NEPC_ext(section,:));
%         hold on
%         polarplot(ELLIPSE_T(section,:),R_ext(section,:));
%         hold off
%         pause();
%         close;
%         
        [ERR_EXT(row),ERR_INT(row)] = get_mean_error(N,npts,NEPC_ext,NEPC_int,R_ext,R_int);
        
    end

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
    
    plot(error_ext(i,:))
    pause();
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


% function [Rnew_ext,Rnew_int] = section_from_PCA(i,ELLIPSE_R_ext,ext_rhoPCAs,ext_rhocoeffs,AVG_RIND_T,numNEPCs)
%     NEPC_ext = zeros(1,size(ext_rhoPCAs,1));
%     
%     for k = 1:numNEPCs
%         % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
%         NEPC_ext = NEPC_ext + ext_rhocoeffs(i,k)*ext_rhoPCAs(:,k)';
%     end
% 
%     Rnew_ext = ELLIPSE_R_ext(i,:) - NEPC_ext;
%     Rnew_int = Rnew_ext - AVG_RIND_T(i);
% 
% end