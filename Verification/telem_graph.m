% analyze OSN telemetry
%% LOAD CSV FILES
% OSN CSV file
telem = readtable('osn_telem.csv');

%% PLOT CSV FILE DATA
figure('Name', 'Telemetry Graph'); 
hold on;

% Plotting Energy vs Time
% Ensure 'timestamp' and 'system_energy' match the headers in your CSV
plot(telem.Var1, telem.Var17, '-r', 'LineWidth', 1.0, 'DisplayName', 'Energy');

% Formatting
grid on; grid minor;
xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold'); % Added essential x-label
ylabel('Energy (J)', 'FontSize', 10, 'FontWeight', 'bold');
title('OSN System Energy Over Time', 'FontSize', 12); % Added title for clarity
legend('show', 'Location', 'southwest');

hold off;
