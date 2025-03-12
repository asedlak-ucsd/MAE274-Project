clear
system = load('fom.mat');
fom = system.sys;


%% Singular Perturbation Model Reduction

sp = SingularPerturbation(fom);

clf;
rom = getrom(sp, 22, false);

sigmaplot(rom);
%%
clf
d = eig(fom.A);
scatter(real(d), imag(d))
hold on;
d = eig(rom.rom.A);
scatter(real(d), imag(d), 'Marker', 'x')
xlabel('$\Re(\lambda)$','interpreter', 'latex')
ylabel('$\Im(\lambda)$','interpreter', 'latex')

%% IRKA Model Reduction 
irka = IRKA(fom);
rom = getrom(irka, 10, 'subspace');


%% Loewner Model Reduction 
n=10;
lh = LoewnerHermite(fom);
sigma = logspace(-1, 4, n);

lh_rom = getrom(lh, sigma, "random", 30);
bodeplot(lh_rom)
