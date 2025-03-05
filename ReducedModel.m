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

        function [e, peak] = error(obj, p)
            % Return error of ROM vs FOM  
            [n, peak] = norm(obj.rom - obj.fom, p);
            e = 100*n / norm(obj.fom, p);
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