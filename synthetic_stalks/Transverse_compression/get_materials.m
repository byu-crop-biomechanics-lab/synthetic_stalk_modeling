function [Erind,Epith] = get_materials(method)
% FILENAME: writespline_V2.m
% AUTHOR: Ryan Larson
% DATE: 5/29/19
%
% PURPOSE: Converts from polar to Cartesian
% 
% 
% INPUTS:
%       method: A string to determine the material selection method.
%       Options include: 
%           'random':   Random material properties
%           'min':      Minimum rind, minimum pith
%           'max':      Maximum rind, maximum pith
%           'minpith':  Random rind, minimum pith
%           'maxpith':  Random rind, maximum pith
%           'minrind':  Minimum rind, random pith
%           'maxrind':  Maximum rind, random pith
%           'avg':      Mean rind, mean pith
% 
% OUTPUTS:
%       Erind: Rind modulus
%       
%       Epith: Pith modulus
%
% NOTES:
%      
% 
% 
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
% Calculate the random material properties from a normal distribution.
    % Bound with 95% confidence interval, calculated from transverse
    % material properties used in another paper.
    Erind_mean = 8.0747e-04; % THESE VALUES ARE IN N/micrometer^2
    Erind_stdev = 3.3517e-04;
    Erind_95 = [6.7414e-04 9.4081e-04];
    Epith_mean = 2.5976e-05;
    Epith_stdev = 1.0303e-05;
    Epith_95 = [2.1878e-05 3.0075e-05];
    ratio_mean = 0.0372;
    ratio_stdev = 0.0180;
    ratio_95 = [0.0300 0.0444];
    
    switch method
        case 'random'
            % Generate Erind from normal distribution
            while 1
                Erind = normrnd(Erind_mean,Erind_stdev);
                if Erind >= Erind_95(1) && Erind <= Erind_95(2)
                    break
                end
            end

            % Generate Epith from normal distribution
            while 1
                Epith = normrnd(Epith_mean,Epith_stdev);
                if Epith >= Epith_95(1) && Epith <= Epith_95(2)
                    break
                end 
            end
    
        case 'min'
            Erind = Erind_95(1);
            Epith = Epith_95(1);
            
        case 'max'
            Erind = Erind_95(2);
            Epith = Epith_95(2);
            
        case 'minpith'
            Erind = Erind_mean;
            Epith = Epith_95(1);
            
        case 'maxpith'
            Erind = Erind_mean;
            Epith = Epith_95(2);
            
        case 'minrind'
            Erind = Erind_95(1);
            Epith = Epith_mean;
            
        case 'maxrind'
            Erind = Erind_95(2);
            Epith = Epith_mean;
    
        case 'avg'
            Erind = Erind_mean;
            Epith = Epith_mean;
    end
    
end