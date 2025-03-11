classdef IRKA
    %IRKA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sys
    end
    
    methods
        function obj = IRKA(sys)
            obj.sys = sys;
        end
        
        function rom = getrom(obj, r, init)

            opts.irka.r = r;
            opts.irka.init = init;
            E = eye(length(obj.sys.A));
            
            [~, Ar, Br, Cr, Dr, ~] = mess_tangential_irka( ...
                E, obj.sys.A, obj.sys.B, obj.sys.C, obj.sys.D, opts);
            
            rom = ReducedModel(obj.sys, ss(Ar, Br, Cr, Dr));
        end
    end
end

