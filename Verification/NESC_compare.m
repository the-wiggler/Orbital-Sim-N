% this file is primarily to plot the telemetry data output by the
% Orbital Sim N (OSN). It has additional capability to plot data from the NESC test
% case CSVs for simulation verification.

%% LOAD CSV FILES
% OSN CSV file
%step001 = readtable('step0-001.csv');
%step0001 = readtable('step0-0001.csv');
%step00001 = readtable('step0-00001.csv');
% Optional NESC file
%NESC_telem = readtable('NESC_test.csv');
telem = readtable('osn_telem.csv');

%% PLOT CSV FILE DATA
figure('Name', 'OSN Telemetry Verification'); 

% --- Subplot 1: X Position ---
subplot(2, 2, 1);
hold on;
plot(telem.timestamp, telem.system_energy, '-r', 'LineWidth', 1.5, 'DisplayName', 'Energy vs Time');
grid on; grid minor;
ylabel('X Position (m)', 'FontSize', 10, 'FontWeight', 'bold');
title('X Position', 'FontSize', 12);
legend('show','Location', 'southwest');
hold off;