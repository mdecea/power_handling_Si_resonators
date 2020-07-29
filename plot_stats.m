function plot_stats(file, T0, Vpp, fig_num, color)

% Plots the data contained in the stats file.

% file --> Filename of the stats file
% T0 --> Operating temeperature
% Vpp --> Peak to peak voltage of the simulated modulation
% fig_num --> The plots will be done starting with figure fig_num
% color --> Color for the lines to plot.

load(file)

OMA = mu_1-mu_0;

figure (fig_num)

subplot(2,2,1)
plot(lamL_m, OMA*1e3, 'Color', color, 'LineWidth', 3)
hold on
xlabel('Wavelength (nm)')
ylabel('OMA (mW)')
set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 30)

subplot(2,2,2)
plot(lamL_m, mu_0*1e3, 'Color', color, 'LineWidth', 3)
hold on
plot(lamL_m, mu_1*1e3, 'Color', color, 'LineWidth', 3)
xlabel('Wavelength (nm)')
ylabel('P_{out} (mW)')
set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 30)

subplot(2,2,3)
plot(lamL_m, T_avg, 'Color', color, 'LineWidth', 3)
hold on
lam_plot = [lamL_m, fliplr(lamL_m)];
coords = [T_max, fliplr(T_min)];
h = fill(lam_plot, coords, 'b', 'LineStyle','none');
h.FaceAlpha=0.3;
h.FaceColor = color;
xlabel('Wavelength (nm)')
ylabel('\DeltaT (K)')
set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 30)

subplot(2,2,4)
plot(lamL_m, N_avg*1e-6, 'Color', color, 'LineWidth', 3)
hold on
lam_plot = [lamL_m, fliplr(lamL_m)];
coords = [N_max*1e-6, fliplr(N_min*1e-6)];
h = fill(lam_plot, coords, 'b', 'LineStyle','none');
h.FaceAlpha=0.3;
h.FaceColor = color;
xlabel('Wavelength (nm)')
ylabel('\DeltaN (cm^{-3})')
set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 30)

set(gcf,'Position',[100 100 1200 (3/4)*1000])

end