T = readtable('orbital_energy.csv');

B1KE = T.Body1_KineticEnergy;
B2KE = T.Body2_KineticEnergy;
B1PE = T.PotentialEnergy;
B2PE = T.PotentialEnergy;

figure;
scatter(T.Distance, B1KE, 1, "red");
hold on;
scatter(T.Distance, B1PE, 1, "blue");
hold on;
scatter(T.Distance, T.TotalOrbitalEnergy, 1, "magenta");
xlabel('Simulation Time');
ylabel('Energy');
legend;
title('Kinetic Energy of Body 1 Over Time');
legend({'Body 1 Kinetic Energy', 'Potential Energy', 'Total Orbital Energy'}, ...
       'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;