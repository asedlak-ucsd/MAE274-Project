classdef SingularPerturbation
   properties
      sys
   end

   methods
       function obj = SingularPerturbation(sys)
            obj.sys = sys;
       end

      function view(obj)
          % Display the participation factors of the system's first 10
          % fastest modes
          [P, d] = participation_factors(obj);
          tiledlayout(2, 5, "TileSpacing", "none");

          for i = 1:10
              nexttile;
              bar(P(i, :));
              ylim([0 1.1]);
              title_text = "$\lambda = " + round(d(i)) + "$";
              title(title_text, 'Position', [length(d) / 2 0.96], 'Interpreter', 'latex');

              set(gca, 'FontSize', 14);
              if i ~= 1 && i ~= 6
                  ax = gca;
                  set(ax, 'YTick', []);
              end

              if i == 1 || i == 6
                ylabel('Participation Factor', 'interpreter','latex', 'FontSize', 16);
                if i == 1
                    ax = gca;
                    set(ax, 'XTick', []);
                end
             end

              if i >= 5
                  xlabel('State', 'Interpreter', 'latex', 'FontSize', 16);
              end
          end
          hold off;
      end

      function rsys = getrom(obj, rstates)
          % Perform a zero-order singular perturbation reduction to 
          % remove the states in rstates from the state-space model

          psys = permute_ss(obj, rstates);

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
   end

   methods (Access = 'protected')
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

      function psys = permute_ss(obj, idx)
          % Permute the order of equations in a state-space model so that all
          % states in idx are the final states
          n = length(obj.sys.A);
          states = 1:n;
          states = [setdiff(states, idx) idx];
          % Generate a permutation matrix P such that P * A will permute the rows 
          % according to the idx array
          P = zeros(n);

          for i = 1:n
              P(i, states(i)) = 1;
          end

          psys = ss2ss(obj.sys, P);
      end
   end
end