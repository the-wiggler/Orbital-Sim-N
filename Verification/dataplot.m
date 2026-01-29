% this file is primarily to plot the telemetry data output by the
% Orbital Sim N (OSN). It has additional capability to plot data from the NESC test
% case CSVs for simulation verification.

%% LOAD CSV FILES
% OSN CSV file
step001 = readtable('step0-001.csv');
step0001 = readtable('step0-0001.csv');
step00001 = readtable('step0-00001.csv');
% Optional NESC file
NESC_telem = readtable('NESC_test.csv');

%% PLOT CSV FILE DATA
figure('Name', 'OSN Telemetry Verification'); 

% --- Subplot 1: X Position ---
subplot(2, 2, 1);
hold on;
plot(step001.timestamp, step001.craft_x, '-r', 'LineWidth', 1.0, 'DisplayName', 'dt=0.001s');
plot(step0001.timestamp, step0001.craft_x, '-g', 'LineWidth', 1.0, 'DisplayName', 'dt=0.0001s');
plot(step00001.timestamp, step00001.craft_x, '-b', 'LineWidth', 1.0, 'DisplayName', 'dt=0.00001s');
plot(NESC_telem.elapsedTime_s, NESC_telem.miPosition_m_X, '-w', 'LineWidth', 1.5, 'DisplayName', 'NESC Ref');
grid on; grid minor;
ylabel('X Position (m)', 'FontSize', 10, 'FontWeight', 'bold');
title('X Position', 'FontSize', 12);
legend('show','Location', 'southwest');
hold off;

% --- Subplot 2: Y Position ---
subplot(2, 2, 2);
hold on;
plot(step001.timestamp, step001.craft_y, '-r', 'LineWidth', 1.0);
plot(step0001.timestamp, step0001.craft_y, '-g', 'LineWidth', 1.0);
plot(step00001.timestamp, step00001.craft_y, '-b', 'LineWidth', 1.0);
plot(NESC_telem.elapsedTime_s, NESC_telem.miPosition_m_Y, '-w', 'LineWidth', 1.5);
grid on; grid minor;
ylabel('Y Position (m)', 'FontSize', 10, 'FontWeight', 'bold');
title('Y Position', 'FontSize', 12);
hold off;

% --- Subplot 3: X Velocity ---
subplot(2, 2, 3);
hold on;
plot(step001.timestamp, step001.craft_vx, '-r', 'LineWidth', 1.0);
plot(step0001.timestamp, step0001.craft_vx, '-g', 'LineWidth', 1.0);
plot(step00001.timestamp, step00001.craft_vx, '-b', 'LineWidth', 1.0);
plot(NESC_telem.elapsedTime_s, NESC_telem.miVelocity_m_s_X, '-w', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('X Velocity (m/s)', 'FontSize', 10, 'FontWeight', 'bold');
title('X Velocity', 'FontSize', 12);
hold off;

% --- Subplot 4: Y Velocity ---
subplot(2, 2, 4);
hold on;
plot(step001.timestamp, step001.craft_vy, '-r', 'LineWidth', 1.0);
plot(step0001.timestamp, step0001.craft_vy, '-g', 'LineWidth', 1.0);
plot(step00001.timestamp, step00001.craft_vy, '-b', 'LineWidth', 1.0);
plot(NESC_telem.elapsedTime_s, NESC_telem.miVelocity_m_s_Y, '-w', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Y Velocity (m/s)', 'FontSize', 10, 'FontWeight', 'bold');
title('Y Velocity', 'FontSize', 12);
hold off;
