classdef SingularPerturbation
   properties
      sys
      epsilon
   end

   methods
       function obj = SingularPerturbation(sys)
            obj.sys = sys;
            % Calculate the minimum (nonzero) absolute value of each row
            % that is, a measure of epsilon in the equation:
            % epsilon*dx = Ax
            obj.epsilon = fast_states(obj);
       end

      function view(obj)
          bar(obj.epsilon)
      end

      function rm = getrom(obj, rstates, order)
          % Perform a zero-order singular perturbation reduction to 
          % remove the n - order fastest states

          psys = permute_ss(obj.sys, rstates);

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
          
          if order == "zero"
              A_inv = inv(A_ff);
              T1 = A_inv * A_fs;
              T2 = A_inv * B_f;
              Ar = A_ss - A_sf * T1;
              Br = B_s - A_sf * T2;
              Cr = C_s - C_f * T1;
              Dr = psys.D - C_f * T2;
    
              rsys = ss(Ar, Br, Cr, Dr);
          end
          if order == "first"
             A_inv = inv(A_ff);
             T1 = A_inv^2 * A_fs;
             T2 = inv(eye(ns) + A_sf*T1);
             T3 = A_sf*A_inv;

             Ar = T2*(A_ss - T3*A_fs);
             Br = T2*(B_s - T3*B_f);

             Cr = C_s - C_f*(T1*Ar + A_inv*A_fs);
             Dr = psys.D - C_f*(T1*Br + A_inv*B_f);

             rsys = ss(Ar, Br, Cr, Dr);
          end

          rm = ReducedModel(obj.sys, rsys);
      end
   end

   methods (Access = 'protected')
       function epsilon = fast_states(obj)

         A = obj.sys.A;
         n = length(A);
         epsilon = zeros(1, n);

         for i=1:n
             % min abs non-zero element per row
             b = A(i,:);
             epsilon(i)=min(abs(b(b ~= 0))); 
         end
      end
   end
end