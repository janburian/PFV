clc
close all
clear all

%% Mereni vzdalenosti pomoci ultrazvuku
distances = [300, 400, 500, 600, 700, 800, 900, 1000]; % [mm]
measured_values_1 = [1.98, 3.03, 4.16, 5.29, 6.42, 7.54, 8.68, 9.76]; % [V]; measured from 300 mm to 1000 mm
measured_values_2 = [1.98, 3.03, 4.16, 5.29, 6.41, 7.5, 8.63, 9.76]; % [V]; measured from 1000 mm to 300 mm (already sorted)

figure;
scatter(distances, measured_values_1, '+', "blue")
hold on
scatter(distances, measured_values_2, '+', "red")
xlabel("d [mm]")
ylabel("U [V]")
title("Mereni vzdalenosti pomoci ultrazvukoveho snimace")


% Prolozeni pomoci primky - Staticka charakteristika
p_1 = polyfit([distances, distances], [measured_values_1, measured_values_2], 1);
x = linspace(300, 1000, 8); % Adapt n for resolution of graph
y = p_1(1) * x + p_1(2);

figure;
scatter(distances, measured_values_1, '+', "blue")
hold on
scatter(distances, measured_values_2, '+', "red")
hold on
plot(x, y, "black"); 
xlabel("d [mm]")
ylabel("U [V]")
title("Staticka charakteristika")
legend("Mereni 1", "Mereni 2", "Aproximacni primka")

% Chyba opakovatelnosti
data_ultra_repetition = readmatrix("./data/Ultra_opak60cm.csv");
data_ultra_repetition_voltage = data_ultra_repetition(:, 4);
%data_ultra_repetition_voltage_avg = sum(data_ultra_repetition(:, 4)) / length(data_ultra_repetition_voltage);
p_2 = polyfit([measured_values_1, measured_values_2], [distances, distances], 1);
x = linspace(300, 1000, 8); % Adapt n for resolution of graph
y = p_2(1) * x + p_2(2);

res = 0;
for i=1:1:length(data_ultra_repetition_voltage)
    res = res + (p_2(1) * data_ultra_repetition_voltage(i) + p_2(2)); 
end
res_avg = res / length(data_ultra_repetition_voltage); 

% p_2 = polyfit(measured_values_1, distances, 1);
% x = linspace(300, 1000, 1000); % Adapt n for resolution of graph
% y = p_2(1) * x + p_2(2);
% 
% figure;
% scatter(distances, measured_values_1, '+', "blue")
% hold on
% scatter(distances, measured_values_2, '+', "red")
% hold on
% plot(x, y, "black"); 
% xlabel("d [mm]")
% ylabel("U [V]")
% title("Chyba opakovatelnosti")

%% Mereni teploty
data_thermometer = readmatrix("./data/teplomer_data_all.csv"); 
data_thermometer_cleaned = data_thermometer(:,[3:5]); % without NaN values [..., temperature, ...]
Ts = 0.1; % perioda vzorkovani

figure
x = linspace(0, length(data_thermometer_cleaned) * Ts, length(data_thermometer_cleaned));
plot(x, data_thermometer_cleaned(:, 2), "-")
hold on 
plot(x, data_thermometer_cleaned(:, 3), "-")
xlabel("Cas [s]")
% ylabel("Teplota [°C]")
title("Vystup referencniho a polovodicoveho snimace")
yline(60, '--', 'Odkryti vodni lazne');
legend("Teplota namerena referencnim snimacem [°C]", "Napeti namerene polovodicovym snimacem U [V]")

% Zavislost teploty na napeti
[max_temp, max_idx] = max(data_thermometer_cleaned(:, 2)); 

p_3 = polyfit(data_thermometer_cleaned(max_idx:end, 2), data_thermometer_cleaned(max_idx:end, 3), 3);
p_4 = polyfit(data_thermometer_cleaned(max_idx:end, 3), data_thermometer_cleaned(max_idx:end, 2), 3);

x = linspace(25, 94, 1000); % Adapt n for resolution of graph
y = p_3(1) * x.^3 + p_3(2) * x.^2 + p_3(3) * x + p_3(4);

figure
plot(data_thermometer_cleaned(max_idx:end, 2), data_thermometer_cleaned(max_idx:end, 3))
hold on
plot(x, y)
xlabel("Teplota namerena referencnim snimacem [°C]")
ylabel("Napeti namerene polovodicovym snimacem U [V]")
title("Zavislost teploty na napeti")
legend("Staticka charakteristika", "Aproximacni polynom")


%% Elektromechanicka soustava modelu pruzne hridele
data_shaft = readmatrix("./data/hridel_ukol_c.csv"); 
data_shaft_cleaned = data_shaft(:,[3:5]); % without NaN and 0 values
Ts = 0.02; % perioda vzorkovani

x = linspace(0, length(data_shaft_cleaned) * Ts, length(data_shaft_cleaned));
figure
plot(x, 10 * sin(Ts * x))
hold on
plot(x, data_shaft_cleaned(:, 2), "-")
hold on 
plot(x, data_shaft_cleaned(:, 3), "-")
%xlim([180 400])
xlabel("Cas [s]")
ylabel("U [V]")

% Staticka charakteristika
data_shaft = readmatrix("./data/hridel_ukol_c.csv"); 
data_shaft_cleaned = data_shaft(:,[3:5]); % without NaN and 0 values
Ts = 0.02; % perioda vzorkovani

[max_shaft, max_idx_shaft] = max(data_shaft_cleaned(:, 2));
[min_shaft, min_idx_shaft] = min(data_shaft_cleaned(:, 2));

x = linspace(-10, 10, length(data_shaft_cleaned(min_idx_shaft:max_idx_shaft, 2)));
figure
plot(x, data_shaft_cleaned(min_idx_shaft:max_idx_shaft, 2));

%% Pruzny pas
data_belt = readmatrix("./data/pruzny_pas_C-a.csv"); 
data_belt_cleaned = data_belt(:,[6:7]); % without NaN and 0 values

figure
plot(data_belt_cleaned(:,1), data_belt_cleaned(:,2))
title("Staticka charakteristika")
xlabel("Laserovy snimac d [mm]")
ylabel("Indukcni snimac U [V]")
hold on
p_5 = polyfit(data_belt_cleaned(:,1), data_belt_cleaned(:,2), 1);
x = linspace(85, 130, 1000); % Adapt n for resolution of graph
y = p_5(1) * x + p_5(2);
plot(x, y)

data_belt_c = readmatrix("./data/pruzny_pas_C-c.csv"); 
data_belt_c_cleaned = data_belt_c(:,[6:7]);



