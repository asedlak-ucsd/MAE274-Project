clear
% Load 5-bus system
system = load('fom2.mat');
fom = system.sys;

A = [-3 -2; 1 0];
B = [1 1; 0 1];
C = [0 1; 1 -1];

sys = ss(A, B, C, zeros(2));

%%
opts.irka.r = 25;
opts.irka.init = 'subspace';

[Er, Ar, Br, Cr, Dr, outinfo] = mess_tangential_irka( ...
    eye(length(fom.A)), fom.A, fom.B, fom.C, fom.D, opts);
rom = ReducedModel(fom, ss(Ar, Br, Cr, Dr))
clf
sigmaplot(rom)
%%

error(rom, Inf)
error(rom, 2)