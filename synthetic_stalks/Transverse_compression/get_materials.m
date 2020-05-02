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
%           'random':           Random material properties
%           'min':              Minimum rind, minimum pith
%           'max':              Maximum rind, maximum pith
%           'minpith':          Mean rind, minimum pith
%           'maxpith':          Mean rind, maximum pith
%           'minrind':          Minimum rind, mean pith
%           'maxrind':          Maximum rind, mean pith
%           'minrind_maxpith':  Minimum rind, maximum pith
%           'maxrind_minpith':  Maximum rind, minimum pith
%           'avg':              Mean rind, mean pith
% 
% OUTPUTS:
%       Erind: Rind modulus
%       
%       Epith: Pith modulus
%
% NOTES:
%      
% -------------------------------------------------------------------------
% SUBROUTINES:
%   N/A
% 
% PSEUDO-CODE:
%   Define mean and standard deviation values for rind stiffness (based on
%   Stubbs 2019 values, in units of N/micrometer^2).
% 
%   Define mean and standard deviation values for pith stiffness (based on
%   Stubbs 2019 values, in units of N/micrometer^2).
% 
%   Depending on the chosen method as defined in the function input,
%   generate the rind and pith material properties. Method descriptions are
%   above in the header.
% 
% -------------------------------------------------------------------------
% VERSION HISTORY:
% V1 - 
% V2 - 
% V3 - 
%
% -------------------------------------------------------------------------
    % Calculate the random material properties from a normal distribution.
    % Bound with 95% confidence interval, calculated from transverse
    % material properties used in another paper (Stubbs 2019).
    Erind_mean = 8.0747e-04; % THESE VALUES ARE IN N/micrometer^2
    Erind_stdev = 3.3517e-04;
    Erind_95 = [6.7414e-04 9.4081e-04];
    Epith_mean = 2.5976e-05;
    Epith_stdev = 1.0303e-05;
    Epith_95 = [2.1878e-05 3.0075e-05];
    
    % Choose which method to use by the input string
    switch method
        % Generate fully random material properties (both rind and pith
        % random) using the mean and standard deviations for rind and pith.
        % Make sure that the values generated are within the 95% confidence
        % interval.
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
        
        % Use lower bound values for both rind and pith
        case 'min'
            Erind = Erind_95(1);
            Epith = Epith_95(1);
            
        % Use upper bound values for both rind and pith
        case 'max'
            Erind = Erind_95(2);
            Epith = Epith_95(2);
            
        % Use lower bound value for rind and upper bound value for pith
        case 'minrind_maxpith'
            Erind = Erind_95(1);
            Epith = Epith_95(2);
            
        % Use upper bound value for rind and lower bound value for pith
        case 'maxrind_minpith'
            Erind = Erind_95(2);
            Epith = Epith_95(1);
            
        % Use lower bound value for pith and mean value for rind
        case 'minpith'
            Erind = Erind_mean;
            Epith = Epith_95(1);
        
        % Use upper bound value for pith and mean value for rind
        case 'maxpith'
            Erind = Erind_mean;
            Epith = Epith_95(2);
        
        % Use lower bound value for rind and mean value for pith
        case 'minrind'
            Erind = Erind_95(1);
            Epith = Epith_mean;
            
        % Use upper bound value for rind and mean value for pith
        case 'maxrind'
            Erind = Erind_95(2);
            Epith = Epith_mean;
        
        % Use mean value for rind and mean value for pith
        case 'avg'
            Erind = Erind_mean;
            Epith = Epith_mean;
    end
    
end