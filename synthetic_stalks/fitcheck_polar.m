% Visually check the fit produced by stalk_cross_fit.m function
close;
fitvals = zeros(990,1);
dmajvals = zeros(990,1);
dminvals = zeros(990,1);
ndepthvals = zeros(990,1);
nwidthvals = zeros(990,1);
nlocvals = zeros(990,1);
aAmpvals = zeros(990,1);
aSymvals = zeros(990,1);

load TRryan.mat rhoDCSR

for i = 1:990
    
    numsection = i; % MAKE SURE THIS IS THE SAME CROSS SECTION THAT WAS OPTIMIZED
    [xopt, fopt, exitflag, output] = stalk_cross_fit_polar(numsection);

    % load polar_cross_sections.mat sections


    % Rdata = sections(numsection,:);
    Rdata = rhoDCSR(1,:,numsection);


    N = length(Rdata);
    theta = linspace(0,2*pi,N);
    % phi = nloc - pi;

    % Get variable values from the output of stalk_cross_fit.m:
    dmaj =      xopt(1);
    dmin =      xopt(2);
    ndepth =    xopt(3);
    nwidth =    xopt(4);
    nloc =      xopt(5);
    aAmp =      xopt(6);
    aSym =      xopt(7);


    % Analysis functions
    asymmetry = aAmp*sin(theta - aSym);
    notch = notch_fn(N,ndepth,nwidth,nloc,theta);
    rsynth = rpts(N,theta,dmaj,dmin,asymmetry,notch);

    % % Plot the original data and the fit data on top of each other
    % polarplot(theta,Rdata);
    % hold on
    % polarplot(theta,rsynth)
    % hold off
    % legend('Original data','Fit model');
    
    fitvals(i) = fopt;
    
    dmajvals(i) =      xopt(1);
    dminvals(i) =      xopt(2);
    ndepthvals(i) =    xopt(3);
    nwidthvals(i) =    xopt(4);
    nlocvals(i) =      xopt(5);
    aAmpvals(i) =      xopt(6);
    aSymvals(i) =      xopt(7);

end

histogram(fitvals,30);
save('RealStalkFit.mat','fitvals','dmajvals','dminvals','ndepthvals','nwidthvals','nlocvals','aAmpvals','aSymvals');



% ------------Stalk geometry functions---------------
function [r] = rpts(N,theta,dmaj,dmin,asymmetry,notch)
    r = zeros(1,N);
    for i = 1:N
        r(i) = (dmaj*dmin/4)/sqrt(((dmin/2)*cos(theta(i)))^2 ...
            + ((dmaj/2)*sin(theta(i)))^2) + asymmetry(i) - notch(i);
    end
end

function [notch] = notch_fn(N,ndepth,nwidth,nloc,theta)
    notch = zeros(1,N);
    for i = 1:N
        notch(i) = ndepth/cosh((10/nwidth)*(theta(i)-nloc))^2;
    end
end


function [fitmetric] = getfitpolar(Rreal,rsynth)
    % Get the overall fit of a synthetic shape against a real shape, in a
    % single value similar to a standard deviation

    N = length(Rreal);

    sum = 0;

    for i = 1:N
        rsq = (rsynth(i) - Rreal(i))^2;
        sum = sum + rsq;
    end

    fitmetric = sqrt(sum/N);

end