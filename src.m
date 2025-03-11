clear
system = load('fom.mat');
fom = system.sys;



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
plot(100*sort(abs(d), 'descend') / max(abs(d)), 'x-')
%ylabel('$H_{\infty}$ Error','interpreter', 'latex')
xlabel('Number of States Removed','interpreter', 'latex')
legend('$H_{\infty}$ Error', '$|\lambda|$','interpreter', 'latex')

%%

%%
irka = IRKA(fom);
rom = getrom(irka, 10, 'subspace');
%%
clf
sigmaplot(rom)
error(rom, Inf)

%%

H = tf(fom);

[p, m] = size(H);
dH_ds = H;   % Preallocate matrix for derivatives

for i = 1:p
    for j = 1:m
        [num, den] = tfdata(H(i,j), 'v');
        dnum = polyder(num);       % Derivative of numerator
        dden = polyder(den);       % Derivative of denominator
        
        % Quotient rule: (f/g)' = (f' * g - f * g') / g^2
        df = tf(dnum, 1);
        dg = tf(dden, 1);
        f = tf(num, 1);
        g = tf(den, 1);
        
        dH_ds(i,j) = (df*g - f*dg) / (g*g);
    end
end

bode(H, dH_ds)

%%
n=7;
lh = LoewnerHermite(fom);
sigma = logspace(-1, 4, n);

lh_rom = getrom(lh, sigma, "random", 30);

bodeplot(lh_rom)

%%
w = logspace(-3, 5, 10000);

clf;
opts = bodeoptions;
opts.PhaseWrapping = 'on';
opts.PhaseWrappingBranch =-180;
bodeplot(a_rom.fom, opts)
hold on;
bodeplot(frd(freqresp(a_rom.rom, w), w), 'r--', opts)
%%
%bodeplot(a_rom)
norm(frd(freqresp(a_rom.rom, w), w) - a_rom.fom, Inf) / norm(a_rom.fom, Inf)

%%
a_rom

%%

% From chat, one way to possibly compute the H2 norm?

% Compute the Frobenius norm (sum of squared singular values) at each frequency
H_frobenius = zeros(size(w));
for i = 1:length(w)
    H_frobenius(i) = norm(lh_rom.rom.ResponseData(:,:,i), 'fro')^2; % Frobenius norm squared
end

% Approximate the positive-frequency H2 norm by numerical integration
H2_pos = sqrt(trapz(w, H_frobenius) / pi);