clear
grid = load('grid.mat');
gfl = load('gfl_src.mat');

sp = SingularPerturbation(gfl.sys);

view(sp)
rsys = getrom(sp, [7 8]);

sigma(gfl.sys, rsys, 'r--')