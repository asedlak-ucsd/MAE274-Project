clear
load('grid.mat')




%% Plot HS and erros
R = lyapchol(sys.A, sys.B);
L = lyapchol(sys.A', sys.C');
s = svd(L'*R);

tiledlayout(1, 2)
nexttile
semilogy(s, '-o');
hold on;
set(gca,'fontsize', 16);
xlabel('Index','interpreter','latex', 'FontSize', 16);
f.Position = [100 100 540 400];
e = 2*cumsum([s; 0],"reverse");
semilogy(e(2:69), '-.diamond');
legend("Hankel Singular Value", "Error Bound", 'interpreter','latex', 'FontSize', 16, 'Location', 'southwest')
hold off;


R = reducespec(sys, "balanced");
R.Options.Algorithm = "absolute";

r=50;
rsys = getrom(R,Order=r);


% Define frequency range
w = logspace(0, 6, 1000);

% Compute singular values for each system
[sv0, ~] = sigma(sys, w);
[sv1, ~] = sigma(rsys, w);
nexttile
semilogx(w, mag2db(sv0),  '-', color=[0 0.4470 0.7410]); hold on;
semilogx(w, mag2db(sv1), 'r--');
set(gca,'fontsize', 16);
hl = legend("FOM (Order 68)", '','','','','','','','','','','','', "ROM (Order "+r+")", 'Location', 'southwest');
set(hl, 'Interpreter','latex', 'FontSize', 16);
ylabel('Singular Values (dB)', 'interpreter', 'latex', FontSize=16)
xlabel('Frequency', 'interpreter', 'latex', 'FontSize', 16)
title(' ')
