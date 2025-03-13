classdef BalancedTruncation
    %IRKA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sys
    end
    
    methods
        function obj = BalancedTruncation(sys)
            obj.sys = sys;
        end
        
        function rom = getrom(obj, r)

            R = reducespec(obj.sys, "balanced");
            R.Options.Algorithm = "relative";           
            rsys = getrom(R,Order=r);
            rom = ReducedModel(obj.sys, rsys);
        end
    end
end

