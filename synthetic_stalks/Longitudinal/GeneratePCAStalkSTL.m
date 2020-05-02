function GeneratePCAStalkSTL(stalknum,npts,nTPCs,nalphaPCs,FileName)
% FILENAME: GeneratePCAStalkSTL.m
% AUTHOR: Ryan Larson
% DATE: 1/31/2020
%
% PURPOSE: 
%   Generate an STL file that is an approximation of a specified stalk used
%   in LongPCAData.mat. Data process is loosely based on
%   GenerateStalkSegmentBetsy_V3.m.
% 
% 
% INPUTS:
%       
%       
% OUTPUTS:
%       - 
%
%
% NOTES: 
%       - 
% 
% 
% VERSION HISTORY:
% V1 - Assumes no meaningful drift (no xbar or ybar components included)
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEV NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/31/2020:
% - There are a lot of different ways the STLs could be generated,
% depending on how many principal components we want to use. The function
% should be changed so it only generates one STL, but it takes the desired
% number of PCs for each variable as inputs. Ideally it should
% automatically name the STLs so it's clear which stalk it matches and how
% many PCs of each variable are used.
% - Need a way to deal with NaN sections of original data. The principal
% components go out past the length of some stalks, so the approximated
% versions need to be cut off somehow to match the original data.
% -STILL NEED TO IMPLEMENT ANGLE CHANGING ALONG THE STALK

% 2/7/2020:
% - Might need to generate rind and pith separately and name them, then
% bring them into Abaqus later.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



load('LongPCAData.mat');

% Set up polar data structures
Theta = zeros(length(keepcols),npts-1); % Rows are slices, columns are points around each cross-section
stalk_ext = zeros(length(keepcols),npts-1);
stalk_int = zeros(length(keepcols),npts-1);

% Verify that stalknum is part of keeprows. Throw an error if it's not a
% keeper.
if ~ismember(stalknum,keeprows)
    error('The chosen stalknum is not part of the longitudinal PCA set. Choose another stalk');
end

% Determine the adjusted index of the stalk of interest
stalkidx = find(keeprows == stalknum);

% Run PCA
[APCs, Acoeffs, ~, ~, ~, ~] = pca(APCA,'Centered',false);
[BPCs, Bcoeffs, ~, ~, ~, ~] = pca(BPCA,'Centered',false);
[TPCs, Tcoeffs, ~, ~, ~, ~] = pca(TPCA,'Centered',false);
[alphaPCs, alphacoeffs, ~, ~, ~, ~] = pca(alphaPCA,'Centered',false);

%% Populate polar data structures
% theta
for i = 1:size(Theta,1)
    theta = linspace(0,2*pi,npts);
    Theta(i,:) = theta(1:end-1);
end

a = Acoeffs(stalkidx,1)*APCs(:,1)'; % Take the first PC for A
b = Bcoeffs(stalkidx,1)*BPCs(:,1)'; % Take the first PC for B

% Make rind thickness vector based on the desired number of PCs included
t = zeros(size(a));
for k = 1:nTPCs
    % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
    t = t + Tcoeffs(stalkidx,k)*TPCs(:,k)';
end

% Make angle vector based on the desired number of PCs included
ang = zeros(size(a));
for k = 1:nalphaPCs
    % Add all NEPCs up to the current NEPC to the ellipse in polar coordinates
    ang = ang + alphacoeffs(stalkidx,k)*alphaPCs(:,k)';
end

% Shift angle vector so node angle orients stalk with minor axis along X
% angshift = mean([ang(18) ang(19)]) + pi/2;
angshift = mean([ang(18) ang(19)]);
for k = 1:length(ang)
    ang(k) = ang(k) - angshift;
end


% Define stalk exterior and interior in polar coordinates (no rotation yet)
for i = 1:size(stalk_ext,1)
    
    theta = Theta(i,:);
    
    r_ext = rpts(npts-1,theta,a(i),b(i));
    
    if all(isnan(r_ext))
        r_int = r_ext;
    else
        r_int = normintV2(r_ext,theta,t(i));
        % Catch cases where the large rind thickness causes errors with the
        % normal offset
        if max(r_int-r_ext) > 0
            r_int = rpts(npts-1,theta,a(i)-2*t(i),b(i)-2*t(i));
        end
    end
    
    stalk_ext(i,:) = r_ext;
    stalk_int(i,:) = r_int;

end

% %% Check cross-sections
% figure(1);
% for i = 1:size(stalk_ext,1)
%     polarplot(Theta(i,:),stalk_ext(i,:),'r');
%     hold on
%     polarplot(Theta(i,:),stalk_int(i,:),'b');
%     hold off
%     pause(0.5);
% end

%% Convert polar data to Cartesian
X_ext = zeros(size(stalk_ext));
Y_ext = zeros(size(stalk_ext));
Z_ext = zeros(size(stalk_ext));
X_int = zeros(size(stalk_int));
Y_int = zeros(size(stalk_int));
Z_int = zeros(size(stalk_int));

slices = mapping(keepcols);

for i = 1:size(stalk_ext,1)
    for j = 1:size(stalk_ext,2)
        X_ext(i,j) = stalk_ext(i,j)*cos(Theta(i,j));
        Y_ext(i,j) = stalk_ext(i,j)*sin(Theta(i,j));
        Z_ext(i,j) = slices(i);
        X_int(i,j) = stalk_int(i,j)*cos(Theta(i,j));
        Y_int(i,j) = stalk_int(i,j)*sin(Theta(i,j));
        Z_int(i,j) = slices(i);
    end
end

%% Rotate and reorder
% Exterior
for i = 1:size(X_ext,1)
    x_ext = X_ext(i,:)';
    y_ext = Y_ext(i,:)';
    
    [~, ~, ~, ~, ~, ~, xp_ext, yp_ext, ~, ~] = reorder_V2_interior(x_ext, y_ext, ang(i), 0, 0);
    
    X_ext(i,:) = xp_ext';
    Y_ext(i,:) = yp_ext';
    
end

% Interior
for i = 1:size(X_int,1)
    x_int = X_int(i,:)';
    y_int = Y_int(i,:)';
    
    [~, ~, ~, ~, ~, ~, xp_int, yp_int, ~, ~] = reorder_V2_interior(x_int, y_int, ang(i), 0, 0);
    
    X_int(i,:) = xp_int';
    Y_int(i,:) = yp_int';
    
end



%% Combine data
flipZ = flip(Z_int);
Zrind = [Z_ext; flipZ];
Zpith = Z_int;

flipX = flip(X_int);
Xrind = [X_ext; flipX];
Xpith = X_int;

flipY = flip(Y_int);
Yrind = [Y_ext; flipY];
Ypith = Y_int;

% Close the rind longitudinally
Xrind = [Xrind; Xrind(1,:)];
Yrind = [Yrind; Yrind(1,:)];
Zrind = [Zrind; Zrind(1,:)];

% Close the rind on each cross-section
Xrind = [Xrind,Xrind(:,1)];
Yrind = [Yrind,Yrind(:,1)];
Zrind = [Zrind,Zrind(:,1)];

% Close the pith longitudinally
Xpith = [zeros(size(Xpith(1,:))); Xpith];
Xpith = [Xpith; Xpith(1,:)];
Ypith = [zeros(size(Ypith(1,:))); Ypith];
Ypith = [Ypith; Ypith(1,:)];
Zpith = [Zpith(1,:); Zpith; Zpith(end,:)];

% Close the pith on each cross-section
Xpith = [Xpith,Xpith(:,1)];
Ypith = [Ypith,Ypith(:,1)];
Zpith = [Zpith,Zpith(:,1)];

FileNameRind = strcat(FileName,'Rind.stl');
FileNamePith = strcat(FileName,'Pith.stl');

surf2stl_V1(FileNameRind,Xrind,Yrind,Zrind);
surf2stl_V1(FileNamePith,Xpith,Ypith,Zpith);


end
