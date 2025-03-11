classdef LoewnerHermite
    %IRKA Summary of this class goes here
    %   Detailed explanation goes here

    properties
        sys
        H
        dH
    end

    methods
        function obj = LoewnerHermite(sys)
            obj.sys = sys;
            obj.H = tf(sys);
            obj.dH = diff(obj);
        end

        function rom = getrom(obj, sigma, directions, n)

            k = length(sigma);
            sigma = [sigma -sigma];

            if directions == "random"
                % Create n ROMs with random directions and return the
                % one with the best H_inf error to FOM.
                err = Inf;
                rom = [];

                for i=1:n
                    % Pick random directions and get ROM
                    r = repmat(rand(2, k), 1, 2);
                    l = repmat(rand(2, k), 1, 2)';
                    % Get Loewner reduced matrices
                    curr_rom = getloewner(obj, sigma, r, l);
                    [curr_err , ~ ] = error(curr_rom, Inf);

                    % Check if ROM is best seen so far and updat if so
                    if curr_err < err
                        err = curr_err;
                        rom = curr_rom;
                    end
                end
            end

            

        end
    end

    methods (Access = 'protected')

        function rom = getloewner(obj, sigma, r, l)
            k = length(sigma);  % Number of columns and rows
            L = zeros(k);       % Loewner matrix
            L_sigma = zeros(k); % Shifted Lowner matrix
            sigma = 1j*sigma;   % Convert to frequencies

            H = freqresp(obj.H, sigma);
            dH = freqresp(obj.dH, sigma);
            W = r;
            V = l;
            
            %%%[1]%%% Form each of the Loewner matrices %%%[1]%%%
            for i = 1:k
                for j = 1:k
                    s_i = sigma(i);
                    s_j = sigma(j);
                    if i == j
                        % Hermite interpolation
                        L(i, i) = l(i,:)*dH(:,:,i)*r(:,i);
                        L_sigma(i, i) = l(i,:)*(H(:,:,i) + s_i*dH(:,:,i))*r(:,i);
                    else
                        % Loewner matrix computation
                        delta = (s_i - s_j);
                        L(i, j) = l(i,:)*(H(:,:,i) - H(:,:,j))*r(:,j) / delta;
                        L_sigma(i, j) = l(i,:)*(s_i*H(:,:,i) - s_j*H(:,:,j))*r(:,j) / delta;
                    end
                end
                W(:, i) = H(:,:,i)*r(:,i);
                V(i, :) = l(i,:)*H(:,:,i);
            end

            %%%[2]%%% Transformation to get real-valued matrices %%%[2]%%%
            T = zeros(k, k);
            for i = 1:k
                if imag(sigma(i)) == 0
                    T(i, i) = 1;
                else
                    [~, j] = min(abs(sigma - conj(sigma(i))));
                    if i < j
                        T(i, i) = 1;
                        T(i, j) = 1;
                        T(j, i) = -1j;
                        T(j, j) = 1j;
                    end
                end
            end

            L = real(T * L * T');
            L_sigma = real(T * L_sigma * T');
            V = real(T * V);
            W = real(W * T');

            %%%[3]%%% Reduce rank if Loewner pencil is singular %%%[2]%%%
            r = rank([L L_sigma]);
            if r < k
                [Y, ~, ~] = svd([L L_sigma]);
                [~, ~, X] = svd([L; L_sigma]);

                Y = Y(:, 1:r);
                X = X(:, 1:r);

                L = Y'*L*X;
                L_sigma = Y'*L_sigma*X;
                V = Y'*V;
                W = W*X;
            end

            rsys = ss(dss(-L_sigma, V, W, 0, -L));
            % Return reduced model object
            rom = ReducedModel(obj.sys, rsys);
        end

        function dH = diff(obj)
            [p, m] = size(obj.H);
            dH = obj.H;

            for i = 1:p
                for j = 1:m
                    [num, den] = tfdata(obj.H(i,j), 'v');
                    dnum = polyder(num);       % Derivative of numerator
                    dden = polyder(den);       % Derivative of denominator

                    % Quotient rule: (f/g)' = (f' * g - f * g') / g^2
                    df = tf(dnum, 1);
                    dg = tf(dden, 1);
                    f = tf(num, 1);
                    g = tf(den, 1);

                    dH(i,j) = (df*g - f*dg) / (g*g);
                end
            end
        end
    end
end
