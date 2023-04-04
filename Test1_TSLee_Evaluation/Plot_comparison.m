clear

model_nobreakage_filename = 'Run1.csv' ;
model_withbreakage_filename = 'EvaluationResults_Run1.csv' ;

model_nobreakage = readtable(model_nobreakage_filename, 'HeaderLines', 1);
model_withbreakage = readtable(model_withbreakage_filename);

figure; hold on; box on;
plot(model_nobreakage.x_m_, model_nobreakage.Hrms_m_, '-b', 'linewidth', 2)
plot(model_withbreakage.x_m_, model_withbreakage.Hrms_m_, '-', 'color', rgb('green'), 'linewidth', 2)
set(gcf,'color','w');
set (gca, 'fontsize', 18)
set(gca, 'linewidth', 2)
xlabel('x (m)')
ylabel('H_{rms} (m)')
legend('without stem breakage', 'with stem breakage')