% this file is primarily to plot the telemetry data output by the
% Orbital Sim N (OSN). It has additional capability to plot data from the NESC test
% case CSVs for simulation verification.


%% LOAD CSV FILES

% OSN CSV file
OSN_telem = readtable('osn_telem.csv');

% Optionoal NESC file (comment out if you just want to see telemetry)
NESC_telem = readtable('NESC_test.csv');

%% PLOT CSV FILE DATA
figure;
plot(OSN_telem.timestamp, OSN_telem.craft_x, '-b', 'LineWidth', 0.5);
hold on;
plot(NESC_telem.elapsedTime_s, NESC_telem.miPosition_m_X, '-r', 'LineWidth', 0.5);
grid on;
xlabel('Time (s)');
ylabel('X Pos (m)');
title('X Pos of Test Craft');