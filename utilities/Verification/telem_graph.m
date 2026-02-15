% OSN CSV file
telem = readtable('osn_telem.csv');

figure('Name', 'Telemetry Analysis'); 

% Create a 2x1 grid (2 rows, 1 column)
t = tiledlayout(2, 1); 

nexttile; 
plot(telem.Var1, telem.Var17, '-r', 'LineWidth', 1.0, 'DisplayName', 'Energy');
title('OSN System Energy Over Time');
ylabel('Energy (J)');
legend
grid on;

nexttile;
plot(telem.Var1, telem.Var3, '-b', 'LineWidth', 1.0, 'DisplayName', 'X Pos');
hold on
plot(telem.Var1, telem.Var4, '-g', 'LineWidth', 1., 'DisplayName', 'Y Pos');
legend;
title('');
ylabel('Value (Units)');
xlabel('Time (s)'); % Usually only need x-label on the bottom plot
grid on;


% Global formatting
title(t, 'OSN Telemetry Overview'); % Overall title for the whole figure