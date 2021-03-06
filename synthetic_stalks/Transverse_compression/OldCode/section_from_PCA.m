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