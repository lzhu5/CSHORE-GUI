clear

model_spatial_filename = 'Run14.csv' ;
model_gage_all_filename = 'var_at_gage7.5.csv' ;

lab_spatial = readtable('lab_measurements_spatial_case14.txt');
lab_gage_all = readtable('lab_measurements_gage_all.txt', 'HeaderLines', 1);

model_spatial = readtable(model_spatial_filename, 'HeaderLines', 1);
model_gage_all = readtable(model_gage_all_filename, 'HeaderLines', 1);

figure; hold on; box on
subplot(121); hold on; box on
plot(model_spatial.x_m_, model_spatial.Hrms_m_, 'k', 'linewidth', 2)
plot(lab_spatial{1, :}, lab_spatial{2, :}, 'or', 'linewidth', 2)
xlim([0, 9.8])
set (gca, 'fontsize', 18)
set(gca, 'linewidth', 2)
set(gcf, 'color', 'w')
xlabel('x (m)')
ylabel('H_{rms} (m)')
legend('model results', 'measurements')
% text(1, 0.0575, '(a)', 'fontsize', 18)

subplot(122); hold on; box on
plot(lab_gage_all.Var2, model_gage_all.Var3, '.r', 'markersize', 18)
fplot(@(x) x, [0, 50], 'k', 'linewidth', 2)
set (gca, 'fontsize', 18)
set(gca, 'linewidth', 2)
set(gcf, 'color', 'w')
xlabel('measured H_{rms} (m)')
ylabel('modeled H_{rms} (m)')
xlim([0, 0.1])
ylim([0, 0.1])
set (gcf, 'Position', [1717         260        1338         485])
axis square
% text(0.08, 0.01, '(b)', 'fontsize', 18)

