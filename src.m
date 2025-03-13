clear
system = load('fom_external.mat');
fom = system.sys;

%% Eigenvalues of the FOM
d = eig(fom.A);
scatter(real(d), imag(d))
xlabel('$\Re(\lambda)$','interpreter', 'latex')
ylabel('$\Im(\lambda)$','interpreter', 'latex')


%%
r = [30 26 22 18 14 10 6 4 2];

% Wrappers for each ROM method
ROM_SP = @(r) getrom(SingularPerturbation(fom), r, true).rom;
ROM_LH = @(r) getrom(LoewnerHermite(fom), logspace(-2, 5, r/2), "random", 25).rom;
ROM_IK = @(r) getrom(IRKA(fom), r, "subspace").rom;
ROM_BT = @(r) getrom(BalancedTruncation(fom), r).rom;

% Compute ROMs each of order R for the four methods
for i=1:length(r)
    sys = ROM_SP(r(i));
    save("sp_"+i, "sys")

    sys = ROM_LH(r(i));
    save("lh_"+i, "sys")

    sys = ROM_IK(r(i));
    save("ik_"+i, "sys")

    sys = ROM_BT(r(i));
    save("bt_"+i, "sys")
end

%% Singular Perturbation 

sp = SingularPerturbation(fom);
view(sp, 10)