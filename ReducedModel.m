classdef ReducedModel
    properties
        fom
        rom
    end

    methods
        function obj = ReducedModel(fom, rom)
            obj.fom = fom;
            obj.rom = rom;
        end

        function e = error(obj, p)
            % Return error of ROM vs FOM  
            e = 100*norm(obj.rom - obj.fom, p) / norm(obj.fom, p);
        end

        function bodeplot(obj)
            opts = bodeoptions;
            opts.PhaseWrapping = 'on';
            opts.PhaseWrappingBranch =-180;
            
            bodeplot(obj.fom, opts);
            hold on;
            bodeplot(obj.rom, 'r--', opts);
            hold off;
        end

        function sigmaplot(obj)
            sigma(obj.fom, obj.rom, 'r--');
        end

        % function for eignvalues

    end

end