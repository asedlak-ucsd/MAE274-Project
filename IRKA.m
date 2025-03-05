classdef IRKA
    %IRKA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sys
        A
        B
        C
    end
    
    methods
        function obj = IRKA(sys)
            obj.sys = sys;
            obj.A = sys.A;
            obj.B = sys.B;
            obj.C = sys.C;
        end
        
        function rom = getrom(obj, sigma, b_hat, c_hat)

            tol = 1e-5;
            converged = 0;
            i = 0;        
            max_iter = 1000;
            while ~converged
                i = i+1;
                Vr = getbasis(obj, sigma, obj.A, obj.B, b_hat);
                Wr = getbasis(obj, sigma, obj.A.', obj.C.', c_hat.');
                
                Ar = Wr.' * obj.A * Vr;

                [X, S, Y] = eig(Ar);
                M = diag(ones(length(sigma), 1) ./ sqrt(diag(Y' * X)));
                Y = Y * M;
                X = X * M;
                b_hat = (Y.' * Wr.'*obj.B).';
                c_hat = (obj.C*Vr * X).';

                prev_sigma = sigma;
                sigma = -1 * diag(S);
                               
                if max(abs(sigma - prev_sigma) ./ abs(prev_sigma)) < tol
                    converged = 1;
                    "Convered in " + i+ " iterations"
                elseif i == max_iter
                    "Failed to converge reached max iterations"
                    converged =1;
                end
            end
            
            rsys = ss(Ar, b_hat.', c_hat.', obj.sys.D);

            rom = ReducedModel(obj.sys, rsys);

        end
    end

    methods (Access = 'protected')
        function U = getbasis(obj, sigma, A, X, y)
            r = length(sigma);
            n = length(A);
            U = zeros(n, r);

            for i=1:r
                U(:, i) = (sigma(i)*eye(n) - A) \ (X*y(:, i));
            end
        end
    end
end

