function [Erind,Epith] = get_materials(method)
% Calculate the random material properties from a normal distribution.
    % Bound with 95% confidence interval, calculated from transverse
    % material properties used in another paper.
    Erind_mean = 8.0747e-04;
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
            
            
    %     % Generate Epith from normal distribution of pith/rind ratios
    %     while 1
    %         ratio = normrnd(ratio_mean,ratio_stdev);
    %         if ratio >= ratio_95(1) && ratio <= ratio_95(2)
    %             break
    %         end
    %     end
    %     Epith = ratio*Erind;

    
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