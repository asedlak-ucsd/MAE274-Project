clear
grid = load('grid.mat');
gfl = load('gfl_src.mat');

sp = SingularPerturbation(grid.sys);

view(sp)

%% ROM for the first GFLI

inputs = [5 6];
outputs = inputs;

rsys = getrom(sp, [1:20 31:68], "first");

opts = bodeoptions;
opts.PhaseWrapping = 'on';
opts.PhaseWrappingBranch =2;

bodeplot(grid.sys(inputs, outputs), opts);
hold on;
bodeplot(rsys(inputs, outputs), 'r--', opts);
hold off;

%% 
sigma(grid.sys(inputs, outputs), rsys(inputs, outputs), 'r--')
100*norm(rsys(inputs, outputs) - grid.sys(inputs, outputs), Inf) / norm(grid.sys(inputs, outputs), Inf)


%% TODO: 
% - Check eigenvalues of ROM
% - Iteratively reducing should not affect ROM, check this