clear
system = load('fom_external.mat');
load('analysis.mat');
load('analysis_ext.mat');
fom = system.sys;

yellow = "#DDAA33";
teal = "#44AA99";
red = "#BB5566";
blue = "#004488";

methods = ["sp", "bt", "ik", "lh"];
colors = [yellow red blue teal];
linestyles = [':', "-.", '--', '-'];
markers = ["square", 'x', 'o', '^'];
names = ["Singular Perturbation", "Balanced Truncation", "Iterative Krylov", "Loewner Hermite"];
%Set default line width and Latex interpeter
set(0, 'DefaultLineLineWidth', 2);
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
   default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
   set(groot, default_name,'latex');
end

clear list_factory system i index_interpreter default_name
%% Singular/Eigenvalue Plot

[G_s, ~] = stabsep(fom);
R = lyapchol(G_s.A, G_s.B)';
L = lyapchol(G_s.A', G_s.C')';
[U, S, V] = svd(L'*R);
s = [Inf; Inf; diag(S)];

d = eig(fom.A);
d = sort(abs(d));

clf
t = tiledlayout(1, 2, "TileSpacing", "tight");
nexttile
semilogy(d, 'Color', yellow, 'LineStyle', ':');
fontsize(gca, 16,"points");
ylabel('$|\lambda(${\boldmath${A}$}$)|$')
xlabel('Index')
yticks([1e0 1e1 1e2 1e3 1e4])

nexttile
semilogy(s, 'Color', red, 'LineStyle', '-.');
fontsize(gca, 16,"points");
ylabel('$\sigma_{\mathcal{H}}$');
xlabel('Index')
yticks([1e-15 1e-10 1e-5 1e0])

%% Error Plot [1]
clf
f = gcf;
f.Renderer = 'painters';

for i=1:4
    data = analysis_external(analysis_external.Method == methods(i), :);
    semilogy(data.Order, data.Hinf, "Color", colors(i), 'Marker', markers(i), ...
        "MarkerSize", 10, "LineStyle", linestyles(i));
    hold on;
end
fontsize(gca, 16,"points");
ylim([1e-11, 2e2])
xlim([0 32])
ylabel("Relative Error (\%)")
xlabel("Model Order");
legend(names, "location",'southwest');
yticks([1e-10 1e-6 1e-2  1e2])

%% Error Plot [2]
clf
f = gcf;
f.Renderer = 'painters';
t = tiledlayout(1, 2);
t.TileSpacing = 'none';

nexttile
for i=1:4
    data = analysis(analysis.Method == methods(i), :);
    semilogy(data.Order, data.Hinf, "Color", colors(i), 'Marker', markers(i), ...
        "MarkerSize", 10, "LineStyle", linestyles(i));
    hold on;
end
fontsize(gca, 16,"points");
ylim([1e-7, 200])
xlim([20 52])
ylabel("Relative Error (\%)")
xlabel("Model Order");
legend(names, "location",'southwest');

nexttile
for i=1:4
    data = analysis(analysis.Method == methods(i), :);
    semilogy(data.Order, data.H2, "Color", colors(i), 'Marker', markers(i), ...
        "MarkerSize", 10, "LineStyle", linestyles(i));
    hold on;
end
fontsize(gca, 16,"points");
ylim([1e-7, 200])
set(gca,'YTickLabel',[]);
xlim([20 52])
xlabel("Model Order")

%% Sigma Plots %%
sys = load("sp_35.mat");
sp_rsys = sys.rsys;
sys = load("bt_35.mat");
bt_rsys = sys.rsys;
sys = load("ik_35.mat");
ik_rsys = sys.rsys;
sys = load("lh_35.mat");
lh_rsys = sys.rsys;

fom_local = load("full_system.mat");
fom_local = fom_local.sys(7:10, 7:11);


clf
f = gcf;
f.Renderer = 'painters';
t = tiledlayout(1, 2);
t.TileSpacing = 'none';
roms = {sp_rsys, bt_rsys, ik_rsys, lh_rsys};
w = logspace(-2, 5, 1000);


nexttile
for i=1:4
    error_sys = roms{i}
    [sv, ~] = sigma(error_sys, w);
    semilogx(w, mag2db(sv), "LineStyle", linestyles(i), "Color", colors(i))
    hold on;
end
[sv, ~] = sigma(fom_local, w);
semilogx(w, mag2db(sv), "LineStyle", "-", "Color", "black", "LineWidth",3)

fontsize(gca, 16,"points");
pad = repmat([""], 1, 3);
legend([names(1) pad names(2) pad names(3) pad names(4) pad "Local FOM"] , "location",'southwest');
set(gca, 'Children', flipud(get(gca, 'Children')));

ylabel("$\sigma$ (dB)");
xlabel("$\omega$ (rad/s)")
ylim([-170 55])
xlim([1, 3e4])
xticks([1e0 1e2 1e4 1e6])


nexttile
for i=1:4
    error_sys = roms{i} - fom_local;
    [sv, ~] = sigma(error_sys, w);
    semilogx(w, mag2db(sv), "LineStyle", linestyles(i), "Color", colors(i))
    hold on;
end
fontsize(gca, 16,"points");
xlim([1, 3e4]);
xticks([1e0 1e2 1e4 1e6])
xlabel("$\omega$ (rad/s)")

ylim([-170 55])
set(gca,'YTickLabel',[]);
%% Eigenvalues in External
clf
f = gcf;
f.Renderer = 'painters';

sys = load("sp_14.mat");
sp_rsys = sys.sys;
sys = load("bt_14.mat");
bt_rsys = sys.sys;
sys = load("ik_14.mat");
ik_rsys = sys.sys;
sys = load("lh_14.mat");
lh_rsys = sys.sys;

roms = {sp_rsys, bt_rsys, ik_rsys, lh_rsys};
box on
for i=1:4
    d = eig(roms{i}.A);
    scatter(real(d), imag(d), 100, "MarkerEdgeColor", colors(i), "Marker", markers(i), "LineWidth", 2)
    box on
    hold on;
end

d = eig(fom.A);
scatter(real(d), imag(d), 100, "MarkerEdgeColor", "black", "Marker", "+", "LineWidth", 2)

legend([names "External FOM"], "Location", "SouthWest")
set(gca, 'Children', flipud(get(gca, 'Children')));
fontsize(gca, 16,"points");
xlabel("$\Re(\lambda)$")
ylabel("$\Im(\lambda)$")
xlim([-3500 100])

% create smaller axes in top left, and plot on it
axes('Position',[.2 .7 .3 .2]);
box on
set(gca,'Box','on');
for i=1:4
    d = eig(roms{i}.A);
    scatter(real(d), imag(d), 100, "MarkerEdgeColor", colors(i), "Marker", markers(i), "LineWidth", 2)
    box on
    hold on;
end
d = eig(fom.A);
scatter(real(d), imag(d), 100, "MarkerEdgeColor", "black", "Marker", "+", "LineWidth", 2)
xlim([-35 5])
ylim([-360 360])

%%
clear
clf

% Set up figure and layout
f = gcf;
f.Renderer = 'painters';
t = tiledlayout(1, 8);
t.TileSpacing = 'none';

% Read data
bar_data = readtable("data_pfs.csv");
cpal = bar_data.ColorGroup;

% First set of bars
ax1 = nexttile([1 3]);
heights = reshape(bar_data.y6(1:30)', 10, 3);
bh1 = bar(heights');

% Apply colors
for i = 1:length(cpal(1:10))
    bh1(i).FaceColor = hex2rgb(cpal{i});
end

% Formatting
fontsize(gca, 16, "points");
set(gca, 'XTickLabel', ["GFLI$_1$" "GFLI$_2$" "GFLI$_3$"]);
ylim([0 1]);
ylabel("$\rho$")

yyaxis right;
yticks([-1 2]);
% Second set of bars
ax2 = nexttile;
heights = reshape(bar_data.y6(31:38), 4, 2);
bh2 = bar(1:4, bar_data.y6(31:34), 'FaceColor', hex2rgb(cpal{31}));
hold on;
bh3 = bar(5:8, bar_data.y6(35:38), 'FaceColor', hex2rgb(cpal{35}));

% Formatting
fontsize(gca, 16, "points");
legend(ax1, [bh1([1 3 5 7 9]) bh2 bh3], {"CC $\frac{1}{s}$", "LCL $i_1$", ...
    "LCL $i_2$", "LCL $v$", "PLL $\frac{1}{s}$", "Shunt $v$", "Line $i$"}, ...
    "Location", "NorthWest", 'NumColumns', 2);

ylim([0 1]); 
ax2.YAxis.Color = 'white';
ax2.YAxis.Visible = 'off';

yticks([-1 2]);
xticks([4.5]);
set(gca, 'XTickLabel', ["Lines"]);


% First set of bars
ax1 = nexttile([1 3]);
heights = reshape(bar_data.y12(1:30)', 10, 3);
bh1 = bar(heights');

% Apply colors
for i = 1:length(cpal(1:10))
    bh1(i).FaceColor = hex2rgb(cpal{i});
end

% Formatting
fontsize(gca, 16, "points");
set(gca, 'XTickLabel', ["GFLI$_1$" "GFLI$_2$" "GFLI$_3$"]);
ylim([0 1]);
set(gca,'YTickLabel',[]);
yyaxis right;
yticks([-1 2]);


% Second set of bars
ax2 = nexttile;
heights = reshape(bar_data.y6(31:38), 4, 2);
bh2 = bar(1:4, bar_data.y12(31:34), 'FaceColor', hex2rgb(cpal{31}));
hold on;
bh3 = bar(5:8, bar_data.y12(35:38), 'FaceColor', hex2rgb(cpal{35}));

% Formatting
fontsize(gca, 16, "points");

ylim([0 1]); 
yyaxis right;
ax = gca;
ax.YAxis(1).Color = 'red';
ax.YAxis(1).Visible = 'off';
ax.YAxis(2).Color = 'k';

yticks([-1 2]);
xticks([4.5]);
set(gca, 'XTickLabel', ["Lines"]);
set ( ax.YAxis(1) , 'visible', 'off' );
print -r600 -dtiff 'test.tif'