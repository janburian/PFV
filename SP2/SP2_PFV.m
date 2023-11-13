clc
close all
clear all

%% Hridel 
% Prechodova charakteristika
data_shaft_step = readmatrix("./data/hridel_step_1.csv"); 
data_shaft_step_cleaned = data_shaft_step(:,[4:7]);
Ts = 0.01; % TODO: spravna perioda vzorkovani? 
x = linspace(0, length(data_shaft_step_cleaned) * Ts, length(data_shaft_step_cleaned));

figure
hold on
plot(x, data_shaft_step_cleaned(:, 3))
plot(x, data_shaft_step_cleaned(:, 1))
plot(x, data_shaft_step_cleaned(:, 2))
xlim([0 100])
title("Prechodova charakteristika")
xlabel("t")
legend("Budici napeti", "Rychlost hridele motoru", "Rychlost hridele setrvacniku")

figure
hold on
plot(x, data_shaft_step_cleaned(:, 3))
plot(x, data_shaft_step_cleaned(:, 1))
plot(x, data_shaft_step_cleaned(:, 2))
ylim([0 4.8])
xlim([19 22])
title("Prechodova charakteristika")
xlabel("t")
legend("Budici napeti", "Rychlost hridele motoru", "Rychlost hridele setrvacniku")


% Impulsni charakteristika
data_shaft_impulse = readmatrix("./data/hridel_imp_1.csv"); 
data_shaft_impulse_cleaned = data_shaft_impulse(:,[4:7]);
Ts = 0.01; % TODO: spravna perioda vzorkovani? 
x = linspace(0, length(data_shaft_impulse_cleaned) * Ts, length(data_shaft_impulse_cleaned));

figure
hold on
plot(x, data_shaft_impulse_cleaned(:, 3))
plot(x, data_shaft_impulse_cleaned(:, 1))
plot(x, data_shaft_impulse_cleaned(:, 2))
ylim([-0.2 0.25])
xlim([20 60])

figure
hold on
plot(x, data_shaft_impulse_cleaned(:, 3))
plot(x, data_shaft_impulse_cleaned(:, 1))
plot(x, data_shaft_impulse_cleaned(:, 2))
ylim([0 0.23])
xlim([23 27])


%% Pruzny pas
% Prechodova charakteristika
data_belt_step = readmatrix("./data/pruzny_pas_step_b2.csv"); 
data_belt_step_cleaned = data_belt_step(:,[4:7]);
Ts = 0.01; % TODO: spravna perioda vzorkovani? 
x = linspace(0, length(data_belt_step_cleaned) * Ts, length(data_belt_step_cleaned));

figure
hold on
plot(x, data_belt_step_cleaned(:, 3))
plot(x, data_belt_step_cleaned(:, 1))
plot(x, data_belt_step_cleaned(:, 2))
xlim([0 100])
title("Prechodova charakteristika")
xlabel("t")
legend("Budici napeti", "Rychlost hridele motoru", "Rychlost hridele setrvacniku")

%% Teplomer





%% Eddy current


