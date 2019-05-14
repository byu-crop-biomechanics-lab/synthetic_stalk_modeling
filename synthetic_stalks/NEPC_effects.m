clear; close;
% load Ellipse_fits_bottom1.mat
load NEPCs_bottom1.mat

%% Individual effects of NEPCs
% Step through the cross sections
for i = 1:length(indices)
    NEPC = zeros(1,size(ext_rhoPCAs,1));
    % Step through the NEPCs
    for j = 1:size(ext_rhocoeffs,1)
        
        % Add the current NEPC to the ellipse in polar coordinates
        NEPC = ext_rhocoeffs(i,j)*ext_rhoPCAs(:,j)'; % reconstruct full scale NEPC for the current cross section
        Rnew_ext = ELLIPSE_R(i,:) - NEPC;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        % Convert from polar to Cartesian
        xy_ext = convert_to_xy(Rnew_ext',ELLIPSE_T(i,:)');
        xy_int = convert_to_xy(Rnew_int',ELLIPSE_T(i,:)');
        xy_original = [ELLIPSE_XY(i,:,1)',ELLIPSE_XY(i,:,2)'];

        % Repeat the last points to close the loop
        xy_ext   = [xy_ext; xy_ext(1,:)];
        xy_int   = [xy_int; xy_int(1,:)];
        xy_original = [xy_original; xy_original(1,:)];
        
        % Plot to check the result
        plot(xy_ext(:,1),xy_ext(:,2));
        hold on
        plot(xy_int(:,1),xy_int(:,2));
        plot(xy_original(:,1),xy_original(:,2),'m');
        axis equal
        titlestring = sprintf('Ellipse %d with PC%d',i,j);
        title(titlestring);
        pause();
        hold off
        
        % Export a .txt file with the new shape data for exterior and
        % interior boundary points
        % Save base points as txt file
        S = size(xy_ext);
        len = max(S);
        outputName_ext = sprintf('Ellipse_wPC%d_ext',j);
        outputName_int = sprintf('Ellipse_wPC%d_int',j);
        % (also export a mat file with the same data for finding the
        % reference points to use in Abaqus)
        matname_ext = strcat(outputName_ext,num2str(i),'.mat');
        matname_int = strcat(outputName_int,num2str(i),'.mat');
        save(matname_ext,'xy_ext');
        save(matname_int,'xy_int');
        writespline(len,xy_ext,i,outputName_ext)
        writespline(len,xy_int,i,outputName_int)
    end
    
    
end


%% Convergence effects of NEPCs

% Step through the cross sections
for i = 1:length(indices)
    NEPC = zeros(1,size(ext_rhoPCAs,1));
    % Step through the NEPCs
    for j = 1:size(ext_rhocoeffs,1)
        
        % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
        NEPC = NEPC + ext_rhocoeffs(i,j)*ext_rhoPCAs(:,j)';
        Rnew_ext = ELLIPSE_R(i,:) - NEPC;
        Rnew_int = Rnew_ext - AVG_RIND_T(i);
        
        % Convert from polar to Cartesian
        xy_ext = convert_to_xy(Rnew_ext',ELLIPSE_T(i,:)');
        xy_int = convert_to_xy(Rnew_int',ELLIPSE_T(i,:)');
        xy_original = [ELLIPSE_XY(i,:,1)',ELLIPSE_XY(i,:,2)'];

        % Repeat the last points to close the loop
        xy_ext   = [xy_ext; xy_ext(1,:)];
        xy_int   = [xy_int; xy_int(1,:)];
        xy_original = [xy_original; xy_original(1,:)];
        
        % Plot to check the result
        plot(xy_ext(:,1),xy_ext(:,2));
        hold on
        plot(xy_int(:,1),xy_int(:,2));
        plot(xy_original(:,1),xy_original(:,2),'m');
        axis equal
        titlestring = sprintf('Ellipse %d with all PCs through PC%d',i,j);
        title(titlestring);
        pause();
        hold off
        
        % Export a .txt file with the new shape data for exterior and
        % interior boundary points
        % Save base points as txt file
        S = size(xy_ext);
        len = max(S);
        outputName_ext = sprintf('Ellipse_wPC1_to_%d_ext',j);
        outputName_int = sprintf('Ellipse_wPC1_to_%d_int',j);
        % (also export a mat file with the same data for finding the
        % reference points to use in Abaqus)
        matname_ext = strcat(outputName_ext,num2str(i),'.mat');
        matname_int = strcat(outputName_int,num2str(i),'.mat');
        save(matname_ext,'xy_ext');
        save(matname_int,'xy_int');
        writespline(len,xy_ext,i,outputName_ext)
        writespline(len,xy_int,i,outputName_int)
    end
    
    
end

close;