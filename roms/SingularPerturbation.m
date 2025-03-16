classdef SingularPerturbation
    properties
        sys
    end

    methods
        function obj = SingularPerturbation(sys)
            obj.sys = sys;
        end

        function view(obj, r)
           % Display the participation factors of the system's first 10
           % fastest modes
           [P, ~] = participation_factors(obj);
           
           % Normalize columns to one, then add the first 38-r rows
           P = P ./ sum(P, 1);
           bar(sum(P(1:(38-r), :), 1));
           ylim([0, 1])
        end

        function rsys = sp_rom(obj, rstates)
            % Perform a zero-order singular perturbation reduction to
            % remove the n - order fastest states

            psys = permute_ss(obj, obj.sys, rstates);

            nf = length(rstates); % ns: number of slow states to include in the ROM
            ns = length(psys.A) - nf;
            A_cells = mat2cell(psys.A, [ns, nf], [ns, nf]);
            B_cells = mat2cell(psys.B, [ns, nf], size(psys.B, 2));
            C_cells = mat2cell(psys.C, size(psys.C, 1), [ns, nf]);

            A_ss = A_cells{1, 1};
            A_sf = A_cells{1, 2};
            A_fs = A_cells{2, 1};
            A_ff = A_cells{2, 2};
            B_s = B_cells{1};
            B_f = B_cells{2};
            C_s = C_cells{1};
            C_f = C_cells{2};

            A_inv = inv(A_ff);
            T1 = A_inv * A_fs;
            T2 = A_inv * B_f;
            Ar = A_ss - A_sf * T1;
            Br = B_s - A_sf * T2;
            Cr = C_s - C_f * T1;
            Dr = psys.D - C_f * T2;

            rsys = ss(Ar, Br, Cr, Dr);
        end

        function rm = getrom(obj, r, feedthrough)
            % Iterative singular perturbation reduction
            rsys = obj.sys;
            
            % Transform system to eigenbasis
            [V, D] = eig(rsys.A);
            tsys = rsys;
            tsys.A = D;
            tsys.B = inv(V)*tsys.B;
            tsys.C = tsys.C*V;

            [~, ids] = sort(diag(abs(D)), "descend");

            % Permute state space model so that states are ordered by
            % the magnitude of eigenvalues
            tsys = permute_ss(obj, tsys, ids);

            % Remove the n fastest states from the current version of the ROM
            n = length(obj.sys.A) - r;
            rsys = sp_rom(SingularPerturbation(tsys), 1:n);

            % Force the resulting ROM to a real valued state space model
            lambda = eig(rsys.A);
            k = length(lambda);
            T = zeros(k, k);

            % Pair up eigenvalues and apply complex transform
            for i = 1:k
                if imag(lambda(i)) == 0
                    T(i, i) = 1;
                else
                    j = i + 1;
                    if conj(lambda(i)) == lambda(j)
                        T(i, i) = 1;
                        T(i, j) = 1;
                        T(j, i) = -1j;
                        T(j, j) = 1j;
                    end
                end
            end

            % Transform the system into real eigenvalue pairs and round 
            % off any remaining small imaginary parts
            try
                rsys = ss2ss(rsys, T);
                rsys.A = real(rsys.A);
                rsys.B = real(rsys.B);
                rsys.C = real(rsys.C);
                rsys.D = real(rsys.D);

                % If the ROM should not have feedthrough set D = 0
                if ~feedthrough
                    rsys.D = rsys.D*0;
                end

                rm = ReducedModel(obj.sys, rsys);
            catch
                % If the reduction removes eigenvalues in an unbalanced
                % manner try increasing the ROM order by one.
                rm = getrom(obj, r+1, feedthrough);
            end
        end
        function [P, d] = participation_factors(obj)
            % Get right eigenvectors and form left eigenvectors using the inverse
            % Ensures the eigenvectors are normalized to one for each mode
            [V, D] = eig(obj.sys.A);
            W = inv(V);
            d = diag(D);
            % Sort the rows of P by each eigenvalue's distance from the origin
            % in the complex plane.
            dP = [abs(d) abs(V .* W')];
            [dP, idx] = sortrows(dP, 'descend');
            P = dP(:, 2:length(d) + 1);
            d = d(idx);
        end

        function psys = permute_ss(obj, sys, idx)
            % Permute the order of equations in a state-space model so that all
            % states in idx are the final states
            n = length(sys.A);
            states = 1:n;
            states = [setdiff(states, idx) idx];
            % Generate a permutation matrix P such that P * A will permute the rows
            % according to the idx array
            P = zeros(n);

            for i = 1:n
                P(i, states(i)) = 1;
            end

            psys = ss2ss(sys, P);
        end
    end
end