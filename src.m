clear
system = load('fom_external.mat');
fom = system.sys;

%% Eigenvalues of the FOM
d = eig(fom.A);
scatter(real(d), imag(d))
xlabel('$\Re(\lambda)$','interpreter', 'latex')
ylabel('$\Im(\lambda)$','interpreter', 'latex')


%%


% Wrappers for each ROM method
ROM_SP = @(r) getrom(SingularPerturbation(fom), r, true);
ROM_LH = @(r) getrom(LoewnerHermite(fom), logspace(-2, 5, r/2), "random", 25);
ROM_IK = @(r) getrom(IRKA(fom), r, "subspace");
ROM_BT = @(r) getrom(BalancedTruncation(fom), r);

%% Compute ROMs each of order R for the four methods
Method = [];
Order = [];
H2 = [];
Hinf = [];
r = [30 26 22 18 14 10 6 4 2];
rom_names = ["bt", "sp", "lh", "ik"];
rom_func = {ROM_BT, ROM_SP, ROM_LH, ROM_IK};

for i=1:4
    for j=1:1%length(r)
        r = 14;
        rom = rom_func{i}(r(j));
        sys = rom.rom;
        save(rom_names(i)+"_"+r, "sys")
        
        Method = [Method; rom_names(i)];
        Order = [Order; length(sys.A)];
        H2 = [H2; error(rom, 2)];
        Hinf = [Hinf; error(rom, Inf)];
    end
end
%%
sigmaplot(rom)
%%
T = table(Method, Order, H2, Hinf)
analysis_external = T;
save("analysis_ext.mat", "analysis_external")
%% Singular Perturbation 

sp = SingularPerturbation(fom);
view(sp, 10)
%%
rom = ROM_LH(10);

%%
sigmaplot(rom)