function GeneratePCAStalkSTL(stalknum,npts,nTPCs,nalphaPCs)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






load('LongPCAData.mat');

% Set up polar data structures
Theta = zeros(length(keepcols),npts); % Rows are slices, columns are points around each cross-section
stalk_ext = zeros(length(keepcols),npts);
stalk_int = zeros(length(keepcols),npts);

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
    Theta(i,:) = linspace(0,2*pi,npts);
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


% Define stalk exterior and interior in polar coordinates (no rotation yet)
for i = 1:size(stalk_ext,1)
    
    theta = Theta(i,:);
    
    r_ext = rpts(npts,theta,a(i),b(i));
    if all(isnan(r_ext))
        r_int = r_ext;
    else
        r_int = normintV2(r_ext,theta,t(i));
    end
    
    stalk_ext(i,:) = r_ext;
    stalk_int(i,:) = r_int;

end


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


%% 
    
end
