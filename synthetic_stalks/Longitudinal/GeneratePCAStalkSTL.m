function GeneratePCAStalkSTL(stalknum,npts)
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

load('LongPCAData.mat');

% Set up polar data structures
Theta = zeros(length(keepcols),npts); % Rows are slices, columns are points around each cross-section
stalk_ext = zeros(length(keepcols),npts);
stalk_int1 = zeros(length(keepcols),npts);
stalk_int2 = zeros(length(keepcols),npts);

% Verify that stalknum is part of keeprows. Throw an error if it's not a
% keeper.
if ~ismember(stalknum,keeprows)
    error('The chosen stalknum is not part of the longitudinal PCA set. Choose another stalk');
end

% Determine the adjusted index of the stalk of interest
stalkidx = find(keeprows == stalknum);

% Run PCA
[APCAs, Acoeffs, ~, ~, ~, ~] = pca(APCA,'Centered',false);
[BPCAs, Bcoeffs, ~, ~, ~, ~] = pca(BPCA,'Centered',false);
[TPCAs, Tcoeffs, ~, ~, ~, ~] = pca(TPCA,'Centered',false);
[alphaPCAs, alphacoeffs, ~, ~, ~, ~] = pca(alphaPCA,'Centered',false);

%% Populate polar data structures
% theta
for i = 1:size(Theta,1)
    Theta(i,:) = linspace(0,2*pi,npts);
end

a = Acoeffs(stalkidx,1)*APCAs(:,1)'; % Take the first PC for A
b = Bcoeffs(stalkidx,1)*BPCAs(:,1)'; % Take the first PC for B
t1 = Tcoeffs(stalkidx,1)*TPCAs(:,1)'; % First PC for T
t2 = Tcoeffs(stalkidx,2)*TPCAs(:,2)'; % Second PC for T
ang1 = alphacoeffs(stalkidx,1)*alphaPCAs(:,1)'; % First PC for alpha
ang2 = alphacoeffs(stalkidx,2)*alphaPCAs(:,2)'; % Second PC for alpha

% Exterior and interior
for i = 1:size(stalk_ext,1)
    
    theta = Theta(i,:);
    
    r_ext = rpts(npts,theta,a(i),b(i));
    if all(isnan(r_ext))
        r_int = r_ext;
    else
        r_int1 = normintV2(r_ext,theta,t1(i));
        r_int2 = normintV2(r_ext,theta,t1(i)+t2(i));
    end
    
    stalk_ext(i,:) = r_ext;
    stalk_int1(i,:) = r_int1;
    stalk_int2(i,:) = r_int2;

end


% THIS BLOCK JUST TAKES THE FITTED A, B, AND T VALUES. IT DOESN'T USE THE
% PRINCIPAL COMPONENTS TO DETERMINE THE VALUES
% % Exterior and interior
% for i = 1:size(stalk_ext,1)
%     dmaj_ext = A(stalknum,i);
%     dmin_ext = B(stalknum,i);
%     rind_t = T(stalknum,i);
%     theta = Theta(i,:);
%     
%     r_ext = rpts(npts,theta,dmaj_ext,dmin_ext);
%     if all(isnan(r_ext))
%         r_int = r_ext;
%     else
%         r_int = normintV2(r_ext,theta,rind_t);
%     end
%     
%     stalk_ext(i,:) = r_ext;
%     stalk_int(i,:) = r_int;
% end

%% Convert polar data to Cartesian


end
