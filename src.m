clear
grid = load('grid.mat');
gfl = load('gfl_src.mat');

sp = SingularPerturbation(grid.sys);

view(sp)

%% ROM for the first GFLI

inputs = [1 2];
outputs = inputs;
sp = SingularPerturbation(grid.sys(inputs, outputs));

sp_rm = getrom(sp, [31:42 47:52 57:64], "zero");



error(sp_rm, Inf)

%% 
bodeplot(sp_rm)


%% TODO: 
% - Check eigenvalues of ROM
% - Iteratively reducing should not affect ROM, check this
% - Perform reduction for the 5-bus system
%   * Remove line from the 5-bus for more radial structure
% Hermite interpolate