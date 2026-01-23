% this file is primarily to plot the telemetry data output by the
% Orbital Sim N (OSN). It has additional capability to plot data from the NESC test
% case CSVs for simulation verification.

%% LOAD CSV FILES
% OSN CSV file
step001 = readtable('step0-001.csv');
step0001 = readtable('step0-0001.csv');
step00001 = readtable('step0-00001.csv');

% Optional NESC file (comment out if you just want to see telemetry)
NESC_telem = readtable('NESC_test.csv');

%% PLOT CSV FILE DATA
figure('Name', 'OSN Telemetry Verification'); 
hold on;

% Simulation Runs (Data)
plot(step001.timestamp, step001.craft_x, ...
    '-r', 'LineWidth', 1.0, 'DisplayName', 'OSN (dt = 0.001s)');

plot(step0001.timestamp, step0001.craft_x, ...
    '-g', 'LineWidth', 1.0, 'DisplayName', 'OSN (dt = 0.0001s)');

plot(step00001.timestamp, step00001.craft_x, ...
    '-b', 'LineWidth', 1.0, 'DisplayName', 'OSN (dt = 0.00001s)');

% NESC Reference
plot(NESC_telem.elapsedTime_s, NESC_telem.miPosition_m_X, ...
    '-w', 'LineWidth', 1.5, 'DisplayName', 'NESC Reference');

grid on;
grid minor;
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('X Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('X-Position Convergence: OSN vs NESC Benchmark', 'FontSize', 14);
lgd = legend('show', 'Location', 'best');
hold off;
