
clear
grid = load('grid.mat');
gfl = load('gfl_src.mat');

function [P, d] = participation_factors(A)
    % Get right eigenvectors and form left eigenvector using inverse
    % This is done to ensure that the eigenvectors are normalized
    % to one for each mode.
    [V, D] = eig(A);
    W = inv(V);
    d = diag(D);
    % Sort the rows of P by each eigenvalues distance from origin
    % in the complex plane.
    dP = [abs(d) abs(V.*W')];
    [dP, idx] = sortrows(dP, 'descend');
    P = dP(:, 2:length(d)+1);
    d = d(idx);
end


function P = permuation_matrix(idx)
    % Generate a permuation matrix P such that P*A will permute the rows 
    % according to the idx array. For example, [2 1 3] will create a matrix
    % that swaps the first and second rows of a 3x3 matrix.

    n = length(idx);
    P = zeros(n);

    for i=1:n
        P(i, idx(i)) = 1;
    end
end


function psys = order_ss(sys, idx)
    % Permuate the order of equations in a state space model
    P = permuation_matrix(idx);
    psys = ss(P'*sys.A*P, P'*sys.B, sys.C*P, sys.D);
end


function rsys = zosp_rom(sys, ns)
    % zero_order_singular_perturbation
    % ns: Number of slow states to include in the ROM
    
    nf = length(sys.A) - ns;
    A_cells = mat2cell(sys.A, [ns, nf], [ns, nf]);
    B_cells = mat2cell(sys.B, [ns, nf], size(sys.B, 2));
    C_cells = mat2cell(sys.C, size(sys.C, 1), [ns, nf]);
    
    A_ss = cell2mat(A_cells(1));
    A_fs = cell2mat(A_cells(2));
    A_sf = cell2mat(A_cells(3));
    A_ff = cell2mat(A_cells(4));
    B_s = cell2mat(B_cells(1));
    B_f = cell2mat(B_cells(2));
    C_s = cell2mat(C_cells(1));
    C_f = cell2mat(C_cells(2));

    A_inv = inv(A_ff);
    T1 = A_inv*A_fs;
    T2 = A_inv*B_f;
    Ar = A_ss - A_sf*T1;
    Br = B_s - A_sf*T2;
    Cr = C_s - C_f*T1;
    Dr = sys.D - C_f*T2;
    
    rsys = ss(Ar, Br, Cr, Dr);
end


%%

sys = order_ss(gfl.sys, [1 2 3 4 5 6 9 10 7 8]);
rsys = zosp_rom(sys, 8);




%%
%%

hold off;
d = eig(rsys.A);
scatter(real(d), imag(d), LineWidth=2)
hold on;
d = eig(sys.A);
scatter(real(d), imag(d), 'x', LineWidth=2, MarkerEdgeColor='r')
hold off;
legend('FOM (Order = 10)', 'ROM (Order = 8)')

%%
bode(gfl.sys, rsys, 'r--');
legend('FOM (Order = 10)', 'ROM (Order = 8)', 'Location','southwest')


%%
[P, d] = participation_factors(rsys.A);
show_pf(P, d)

%%










%% BT Code
R = lyapchol(A, B)';
L = lyapchol(A', C')';
[U, S, V] = svd((L')*R);
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);
T = R*Vr*Sr^(-0.5);
T_inv = Sr^(-0.5)*(Ur')*(L');




















x_states = ["IntCurrController_d"
            "IntCurrController_q"
            "LCLcurrent1_d"
            "LCLcurrent1_q"
            "LCLcurrent2_d"
            "LCLcurrent2_q"
            "LCLvoltcapacitor_d"
            "LCLvoltcapacitor_q"
            "IntPLL"
            "PhasePLL"];

function show_pf(P, d)
    tiledlayout(2, 4, "TileSpacing","none");
    
    for i = 1:10
        nexttile
        
        bar(P(i, :));
        ylim([0 1.09]);
        title_text = "$\lambda = "+round(d(i))+"$";
        title(title_text, 'Position',[5 0.96], 'interpreter','latex');
        
        set(gca,'fontsize', 14);
        if i ~= 1
            if i~=5
                ax = gca;
                set(ax, 'YTick', []);
            end
        end

        if i == 1 | i == 5
            ylabel('Participation Factor', 'interpreter','latex', 'FontSize', 16);
            if i == 1
                ax = gca;
                set(ax, 'XTick', []);
            end
        end

        if i>=5
            xlabel('State', 'interpreter','latex', 'FontSize', 16);
        end
    end
    hold on;
end


[P, d] = participation_factors(gfl.sys.A);
show_pf(P, d)


%%
%%
e_fom = eig(sys.A);
scatter(real(e_fom), imag(e_fom));
set(gca,'fontsize', 14);
xlabel('$\Re(\lambda)$', 'interpreter','latex', 'FontSize', 16);
ylabel('$\Im(\lambda)$', 'interpreter','latex', 'FontSize', 16, 'Rotation', 0);

%%
function rsys = reduce_states(A, B, C, D, states)
    
    [Ar, Br, Cr, Dr] = row_reduce(A, B, C, D, states(1));
    for k=2:length(states)
        i = states(k);
        [Ar, Br, Cr, Dr] = row_reduce(Ar, Br, Cr, Dr, i-k);
    end

    rsys = ss(Ar, Br, Cr, Dr);
end


function [Ar, Br, Cr, Dr] = row_reduce(A, B, C, D, i)
    % Eliminate state i from all rows of A by assuming dx_i = 0 that is
    % 0 = A(i, :)*x + B(i, :)*u
    
    % The amount to scale the i-th row of A by to eliminate x_i
    alpha = A(:, i) ./ A(i, i);
    beta = C(:, i) ./ A(i, i);

    % Eliminate x_i from the j-th row of A and B
    Ar = A - alpha .* A(i, :);
    Br = B - alpha .* B(i, :);

    % Eliminate x_i from C and D
    Cr = C - beta .* A(i, :);
    Dr = D - beta .* B(i, :);

    Ar(:, i) = [];
    Ar(i, :) = [];
    Br(i, :) = [];
    Cr(:, i) = [];
end

rsys = reduce_states(gfl.sys.A, gfl.sys.B, gfl.sys.C, gfl.sys.D, [7]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [3]);

%%
tiledlayout(1, 2, "TileSpacing","compact");
% Define frequency range
w = logspace(0, 6, 1000);

% Compute singular values for each system
[sv0, ~] = sigma(gfl.sys, w);
[sv1, ~] = sigma(rsys, w);
nexttile
semilogx(w, mag2db(sv0),  '-', color=[0 0.4470 0.7410]); hold on;
semilogx(w, mag2db(sv1), 'r--');
set(gca,'fontsize', 16);
hl = legend("FOM (Order 10)",'',"ROM (Order 8)", "Location", 'southwest');
set(hl, 'Interpreter','latex', 'FontSize', 16);
ylabel('Singular Values (dB)', 'interpreter', 'latex', FontSize=16)
xlabel('Frequency', 'interpreter', 'latex', FontSize=16)
title(' ')


% Dropping states with greates contribution to the most neg eigenvalue
rsys = reduce_states(grid.sys.A, grid.sys.B, grid.sys.C, grid.sys.D, [55]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [57 58]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [53 54]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [33]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [34]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [56]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [35]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [34 33]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [13]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [31]);

% Define frequency range
w = logspace(0, 6, 1000);

% Compute singular values for each system
[sv0, ~] = sigma(grid.sys, w);
[sv1, ~] = sigma(rsys, w);
nexttile
semilogx(w, mag2db(sv0),  '-', color=[0 0.4470 0.7410]); hold on;
semilogx(w, mag2db(sv1), 'r--');
set(gca,'fontsize', 16);
hl = legend("FOM (Order 68)",'','','','','','','','','','','','',"ROM (Order 55)");
set(hl, 'Interpreter','latex', 'FontSize', 16, "Location", 'southwest');
ylabel('Singular Values (dB)', 'interpreter', 'latex', FontSize=16)
xlabel('Frequency', 'interpreter', 'latex', 'FontSize', 16)
title(' ')
%%
[P, d] = participation_factors(sys.A);
show_pf(P, d)

%%
rsys = reduce_states(sys.A, sys.B, sys.C, sys.D, [23 31])% 31 32]);
sigma(sys, rsys, 'r--')
%%
rsys = reduce_states(sys.A, sys.B, sys.C, sys.D, [58, 59]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [56, 57]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [54, 55]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [53, 52]);
rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [50, 51]);
sigma(sys, rsys, 'r--')

%%
rsys = reduce_states(sys.A, sys.B, sys.C, sys.D, [17, 18]);
for i=1:6
    k=2;
    rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [17-k*i, 18-k*i]);
end

%for i=1:20
%    rsys = reduce_states(rsys.A, rsys.B, rsys.C, rsys.D, [63-i]);
%end

sigma(sys, rsys, 'r--')
%%
P = participation_factors(rsys.A);
show_pf(P)

%%
e_fom = eig(sys.A);
e_rom =eig(rsys.A);
scatter(real(e_fom), imag(e_fom)) 
hold on;
scatter(real(e_rom), imag(e_rom), 'red')
hold off;
%%
A = magic(3);
B = A(:, 2:3);
C = magic(4).*2;
C = C(:, 2:4);
D = zeros(4);
D = D(:, 1:2);


[Ar, Br, Cr, Dr] = row_reduce(A, B, C, D, 1)


%%
function rsys = eliminate_states(sys, states)
    nx = size(sys.B, 1);
    nu = size(sys.B, 2);
    ny = size(sys.C, 1);

    E = eye(nx);
    for k=1:length(states)
        i = states(k);
        E(i, i) = 0;
    end
    
    x = sym("x", [nx, 1]);
    dx = sym("dx", [nx, 1]);
    u = sym("u", [nu, 1]);
    y = sym("y", [ny, 1]);


    eqs1 = E * dx == sys.A * x + sys.B * u;
    eqs2 = y == sys.C * x + sys.D * u;
    
    for k=1:length(states)
        i = states(k);
        % Eliminate x_i frim dx = Ax + Bu
        eqs1 = eliminate(eqs1, x(i));
        % Solve for x_i and plug it into y = Cx + Du
        xi_solved = solve(eqs1(i), x(i));
        eqs2 = subs(eqs2, x(i), xi_solved);
    end

    % Convert to matrix form
    Er = equationsToMatrix(eqs1, dx);
    Ar = equationsToMatrix(eqs1, x);
    Br = equationsToMatrix(eqs1, u);
    Cr = equationsToMatrix(eqs2, x);
    Dr = equationsToMatrix(eqs2, u);
    
    for k=1:length(states)
        i = states(k);
        Er(:, i) = [];
        Ar(:, i) = [];
        Cr(:, i) = [];
    end

    dsys = dss(Ar, Br, Cr, Dr, Er);
    rsys = dss2ss(dsys);
end


rsys = eliminate_states(sys, [7, 8]);

%%
A = [-1.5,-2;1,0];
B = [0.5;0];
C = [0,1];
D = 0;
sys = ss(A,B,C,D);
rsys = eliminate_states(sys, []);
%%
red = eliminate(eqs, x2)
equationsToMatrix(red, u)
%%
clear
syms x1 x2 x3 u1 u2 u3


% Define state and derivative vectors
x = [x1; x2; x3];
u = [u1; u2; u3];
dx = diff(x);

% Define E, A, and B matrices
E = [1  0  0; 
     0  0  0;
     0  0  1]; % Descriptor matrix with an algebraic constraint

A = [1 1 0; 
     4 1 1;
     3 0 2];

B = eye(3);

% Descriptor system equation: E * dx = A * x + B * u
eqs = 0 == A * x + B * u;

[newEqs, newVars] = reduceRedundancies(eqs, x);
%%
% Step 2: Identify and solve for algebraic constraint
% Find equation(s) with no derivatives (purely algebraic)
alg_eq = newEqs(2); % Example: picking the second equation if it's algebraic
solve_for_x2 = solve(alg_eq, x2); % Solve for x2 in terms of other variables

% Step 3: Substitute into remaining equations
newEqs = subs(newEqs, x2, solve_for_x2);

% Extract matrices for state-space form
[Ar, b] = equationsToMatrix(newEqs, newVars);
Br = equationsToMatrix(b, u);
%%

