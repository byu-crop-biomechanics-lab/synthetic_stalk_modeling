function plotfitstalk_polar(N,dmaj,dmin,ndepth,nwidth,nloc,aAmp,aSym)
% plotfitstalk_polar.m: Quick plotting function for plotting the polar
% version of the fitting equations (as of 4/19/2019) to test the effects of
% different parameter values.
theta = linspace(0,2*pi,N);

asymmetry = aAmp*sin(theta - aSym);
notch = notch_fn(N,ndepth,nwidth,nloc,theta);
r = rpts(N,theta,dmaj,dmin,asymmetry,notch);

polarplot(theta,r);

end