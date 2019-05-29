function getfitpts(cross_section,csvname)
%getfitpts.m: Create a CSV file of the points that define a chosen cross
%section that was fit from real data


% Load parameter values
% clear;
load RealStalkFit.mat aAmpvals aSymvals dmajvals dminvals ndepthvals nlocvals nwidthvals

% nsections = 990;
Npts = 90;

theta = linspace(0,2*pi,Npts);

% Set parameter values for current cross section
dmaj =  dmajvals(cross_section);
dmin =  dminvals(cross_section);
ndepth =  ndepthvals(cross_section);
nwidth = nwidthvals(cross_section);
nloc = nlocvals(cross_section);
aAmp = aAmpvals(cross_section);
aSym = aSymvals(cross_section);

% Calculate cross section shape
asymmetry = aAmp*sin(theta - aSym);
notch = notch_fn(Npts,ndepth,nwidth,nloc,theta);
r = rpts(Npts,theta,dmaj,dmin,asymmetry,notch);

% Convert from polar to Cartesian coordinates (must include Z coordinates
% for the Fusion 360 script to convert properly)
X = zeros(Npts,1);
Y = zeros(Npts,1);
Z = zeros(Npts,1);

for i = 1:Npts
    X(i) = r(i)*cos(theta(i));
    Y(i) = r(i)*sin(theta(i));
end

data = [X,Y,Z];

csvwrite(csvname,data);

end

