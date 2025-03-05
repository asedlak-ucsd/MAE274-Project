clear
% Load 5-bus system
system = load('fom.mat');
fom = system.sys;


%% BT Reductions

% Percent reduction in states: 84.21%, 78.95%, 73.68%, 50%
r = [5, 11, 19, 23, 28];

for i=1:length(r) 
    R = reducespec(fom, "balanced");
    R.Options.Algorithm = "relative";
    sys = getrom(R, Order=r(i));
    save("rom_bt"+r(i), "sys");
end

%% SP Reductions

% Greedy algorithm to find the best states to remove
sp = SingularPerturbation(fom);

r = 35;
s = [];
errors = [];
for i=1:r
    states = setdiff(1:38, s);
    min_error = Inf;
    min_error_state = 0;
    for j=1:r-i
        try
            sp_rom = getrom(sp, [s states(j)], "first");
            curr_error = error(sp_rom, Inf);
        catch exception
            continue;
        end        

        if curr_error < min_error
            min_error = curr_error;
            min_error_state = states(j);
        end
    end
    s = [s min_error_state];
    errors = [errors min_error];
end
%%
d = eig(fom.A);
scatter(real(d), imag(d))
xlabel('$\Re(\lambda)$','interpreter', 'latex')
ylabel('$\Im(\lambda)$','interpreter', 'latex')
%%
clf
plot(errors, 'o-')
hold on;
plot(100*sort(abs(real(d))) / max(abs(real(d))), 'x-')
%ylabel('$H_{\infty}$ Error','interpreter', 'latex')
xlabel('Number of States Removed','interpreter', 'latex')
legend('$H_{\infty}$ Error', '$|\Re(\lambda)|$','interpreter', 'latex')

%%
clf
sp = SingularPerturbation(fom);
rom = getrom(sp, s(1:14), "first");
d = eig(fom.A);
scatter(real(d), imag(d))
xlabel('$\Re(\lambda)$','interpreter', 'latex')
ylabel('$\Im(\lambda)$','interpreter', 'latex')
hold on;
d = eig(rom.rom.A);
scatter(real(d), imag(d), 'red', 'x')
xlabel('$\Re(\lambda)$','interpreter', 'latex')
ylabel('$\Im(\lambda)$','interpreter', 'latex')
%%
sp = SingularPerturbation(fom);
rom = getrom(sp, s(1:14), "first");
sigmaplot(rom);
led = legend('FOM (Order=38)', 'First-Order SP (Order=24)');
set(led, 'interpreter', 'latex');
%%
sp = SingularPerturbation(fom);
rom = getrom(sp, s(1:20), "first");
sigmaplot(rom);
led = legend('FOM (Order=38)', 'First-Order SP (Order=20)');
set(led, 'interpreter', 'latex');
%%
sp = SingularPerturbation(fom);
rom = getrom(sp, s(1:28), "first");
sigmaplot(rom);
led = legend('FOM (Order=38)', 'First-Order SP (Order=10)');
set(led, 'interpreter', 'latex');
%%



%%



r = [27, 33];

for i=1:length(r)
    sp_rom = getrom(sp, s(1:r(i)), "first");
    sys = sp_rom.rom;
    save("rom_fsp"+(38-r(i)), "sys");
end



%% TODO: 
% - Hermite interpolate
% - What frequencies are relevant in dq-frame


%%
[D, V] = eig(fom.A);
U = zeros(38);

for i=1:38
    if i == 1
        U(:, i) = real(V(:, i));
    elseif real(V(:, i)) == real(V(:, i-1)) 
        U(:, i) = imag(V(:, i));
    else
        U(:, i) = real(V(:, i));
    end
end


out = ss2ss(fom, inv(U));

%out.A = real(out.A)

%%
sigma(fom, out)


%%
sp = SingularPerturbation(out);
r = 35;
s = [];
errors = [];
for i=1:r
    states = setdiff(1:38, s);
    min_error = Inf;
    min_error_state = 0;
    for j=1:r-i
        try
            sp_rom = getrom(sp, [s states(j)], "first");
            curr_error = error(sp_rom, Inf);
        catch exception
            continue;
        end        

        if curr_error < min_error
            min_error = curr_error;
            min_error_state = states(j);
        end
    end
    s = [s min_error_state];
    errors = [errors min_error];
end
%%

rom = getrom(sp, s(1:14), 'first');

error(rom, Inf)
sigmaplot(rom)

%%
d = eig(rom.rom.A);

scatter(real(d), imag(d))%

%plot(abs(d))