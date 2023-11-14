clc
close all
clear all

%% Hridel 
% Prechodova charakteristika
data_shaft_step = readmatrix("./data/hridel_step_1.csv"); 
data_shaft_step_cleaned = data_shaft_step(:,[4:7]);
Ts = 0.01; % TODO: spravna perioda vzorkovani? 
x = linspace(0, length(data_shaft_step_cleaned) * Ts, length(data_shaft_step_cleaned));
U_max = 5;
normalization_step = 1 / U_max;

data_shaft_step_cleaned_normalized = data_shaft_step_cleaned * normalization_step;

figure
hold on
plot(x, data_shaft_step_cleaned_normalized(:, 3))
plot(x, data_shaft_step_cleaned_normalized(:, 2))
xlim([0 78])
title("Prechodova charakteristika")
xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure
plot(x, data_shaft_step_cleaned_normalized(:, 1))
xlim([0 78])
title("Prechodova charakteristika hridele motoru")
xlabel("t")
legend("Rychlost hridele motoru")

figure
hold on
plot(x, data_shaft_step_cleaned_normalized(:, 3))
plot(x, data_shaft_step_cleaned_normalized(:, 2))
ylim([0 1])
xlim([19 22])
title("Prechodova charakteristika")
xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")

% Impulsni charakteristika
data_shaft_impulse = readmatrix("./data/hridel_imp_1.csv"); 
data_shaft_impulse_cleaned = data_shaft_impulse(:,[4:7]);
Ts = 0.01; % TODO: spravna perioda vzorkovani? 
U_max = 5;
T_max = 10;
normalization_impulse = 1 / (U_max + T_max); % TODO: normalizace
x = linspace(0, length(data_shaft_impulse_cleaned) * Ts, length(data_shaft_impulse_cleaned));

data_shaft_impulse_cleaned_normalized = data_shaft_impulse_cleaned * normalization_impulse;

figure
hold on
plot(x, data_shaft_impulse_cleaned_normalized(:, 3))
%plot(x, data_shaft_impulse_cleaned_normalized(:, 1))
plot(x, data_shaft_impulse_cleaned_normalized(:, 2))
title("Impulsni charakteristika")
% ylim([-0.2 0.25])
% xlim([20 60])

% figure
% hold on
% plot(x, data_shaft_impulse_cleaned(:, 3))
% plot(x, data_shaft_impulse_cleaned(:, 1))
% plot(x, data_shaft_impulse_cleaned(:, 2))
% ylim([0 0.23])
% xlim([23 27])

% Frekvencni charakteristika TODO: vyextrahovat frekvence 
data_shaft_freq = readmatrix("./data/hridel_frek_1.csv"); 
data_shaft_freq_cleaned = data_shaft_freq(:,[4:7]);

x = linspace(0, length(data_shaft_freq_cleaned(2000:3000)) * Ts, length(data_shaft_freq_cleaned(2000:3000)));

figure 
hold on
plot(x, data_shaft_freq_cleaned(2000:3000, 3)) % Budici napeti
plot(x, data_shaft_freq_cleaned(2000:3000, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika 1")
%xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure 
hold on
plot(x, data_shaft_freq_cleaned(6000:7000, 3)) % Budici napeti
plot(x, data_shaft_freq_cleaned(6000:7000, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika 2")
%xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure 
hold on
plot(x, data_shaft_freq_cleaned(11000:12000, 3)) % Budici napeti
plot(x, data_shaft_freq_cleaned(11000:12000, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika 3")
%xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")


%% Pruzny pas
% Prechodova charakteristika
data_laser_belt_step = readmatrix("./data/pruzny_pas_step_b1.csv");
data_induction_belt_step = readmatrix("./data/pruzny_pas_step_b2.csv"); 

data_belt_step_cleaned = data_induction_belt_step(:,[4:7]);

Ts = 0.01; % TODO: spravna perioda vzorkovani? 
x = linspace(0, length(data_belt_step_cleaned) * Ts, length(data_belt_step_cleaned));
U_max = 5;
normalization_step = 1 / U_max;

data_belt_step_cleaned_normalized = data_belt_step_cleaned * normalization_step;

figure
hold on
plot(x, data_belt_step_cleaned_normalized(:, 4)) 
%plot(x, data_belt_step_cleaned_normalized(:, 1))
%plot(x, data_belt_step_cleaned_normalized(:, 2))
%xlim([0 100])
title("Prechodova charakteristika")
xlabel("t")
legend("Budici napeti", "Rychlost hridele motoru", "Rychlost hridele setrvacniku")

%% Teplomer
data_thermometer = readmatrix("./data/teplomer_data_all.csv"); 
data_thermometer_cleaned = data_thermometer(:,[3:5]); % without NaN values [..., temperature, ...]
Ts = 0.1; % perioda vzorkovani
U_max = 10;

data_thermometer_25_to_90_temperature = data_thermometer_cleaned(1:max_idx, 2);
data_thermometer_90_to_25_temperature = data_thermometer_cleaned(max_idx:end, 2);
data_thermometer_25_to_90_voltage = data_thermometer_cleaned(1:max_idx, 3);

[max_temp, max_idx] = max(data_thermometer_25_to_90_temperature); 
[min_temp, min_idx] = min(data_thermometer_25_to_90_temperature); 

normalization_step_temp = 1 / (100 - 0);
normalization_step_voltage = 1 / U_max;

x = linspace(0, length(data_thermometer_25_to_90_temperature) * Ts, length(data_thermometer_25_to_90_temperature));

data_thermometer_25_to_90_temperature_normalized = data_thermometer_25_to_90_temperature * normalization_step_temp;
data_thermometer_25_to_90_voltage_normalized = data_thermometer_25_to_90_voltage * normalization_step_voltage;

figure
hold on
plot(x, data_thermometer_25_to_90_temperature_normalized)
plot(x, data_thermometer_25_to_90_voltage_normalized)
xlim([345 380])
% xlabel("Teplota namerena referencnim snimacem [°C]")
% ylabel("Napeti namerene polovodicovym snimacem U [V]")
title("Prechodova charakteristika")
% legend("Staticka charakteristika", "Aproximacni polynom")

%% Eddy current
% Prechodova charakteristika
data_eddy_step = readmatrix("./data/Eddy_step.csv"); 
data_eddy_step_cleaned = data_eddy_step(:,[4:7]);
Ts = 0.01; % TODO: spravna perioda vzorkovani? 
x = linspace(0, length(data_eddy_step_cleaned) * Ts, length(data_eddy_step_cleaned));
U_max = 5;
normalization_step = 1 / U_max;

figure
hold on
plot(x, data_eddy_step_cleaned(:, 3))
%plot(x, data_eddy_step_cleaned(:, 2))
plot(x, data_eddy_step_cleaned(:, 1))
xlim([0 13])
title("Prechodova charakteristika")
xlabel("t")

% Frevencni charakteristika
data_eddy_freq = readmatrix("./data/Eddy_frek.csv"); 
data_eddy_freq_cleaned = data_eddy_freq(:,[4:7]);

x = linspace(0, length(data_eddy_freq_cleaned(1:400)) * Ts, length(data_eddy_freq_cleaned(1:400)));

figure
hold on
plot(x, data_eddy_freq_cleaned(1:400, 3))
%plot(x, data_eddy_freq_cleaned(:, 2))
plot(x, data_eddy_freq_cleaned(1:400, 1))
title("Frekvencni charakteristika")
xlabel("t")

