clc
close all
clear all

%% Hridel 
% Prechodova charakteristika
data_shaft_step = readmatrix("./data/hridel_step_1.csv"); 
data_shaft_step_cleaned = data_shaft_step(:,[4:7]);
Ts = 0.005;
x = linspace(0, length(data_shaft_step_cleaned) * Ts, length(data_shaft_step_cleaned));
U_max = 5;
normalization_step = 1 / U_max;

data_shaft_step_cleaned_normalized = data_shaft_step_cleaned * normalization_step;

figure
hold on
plot(x, data_shaft_step_cleaned_normalized(:, 3))
plot(x, data_shaft_step_cleaned_normalized(:, 2))
xlim([0 50])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure
plot(x, data_shaft_step_cleaned_normalized(:, 1))
xlim([0 50])
title("Prechodova charakteristika hridele motoru")
xlabel("t [s]")
ylabel("y(t)")
legend("Rychlost hridele motoru")
%%
ipt = findchangepts(data_shaft_step_cleaned_normalized(1:2100, 2));
r=-data_shaft_step_cleaned_normalized(ipt+1, 2)+data_shaft_step_cleaned_normalized(ipt-1, 2);
p1=-[10*(x(ipt-1)-x(ipt+1))-x(ipt);10*r+data_shaft_step_cleaned_normalized(ipt)];
p2= [10*(+x(ipt-1)-x(ipt+1))+x(ipt);10*r+data_shaft_step_cleaned_normalized(ipt)];

figure
hold on
plot(x, data_shaft_step_cleaned_normalized(:, 3))
plot(x, data_shaft_step_cleaned_normalized(:, 2))
plot(x(ipt),data_shaft_step_cleaned_normalized(ipt,2),'ko');
plot([p1(1) p2(1)],[p1(2) p2(2)],'k--')
ylim([0 1])
xlim([19 22])
title("Prechodova charakteristika")
xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")

%%
% Impulsni charakteristika
data_shaft_impulse = readmatrix("./data/hridel_imp_1.csv"); 
data_shaft_impulse_cleaned = data_shaft_impulse(:,[4:7]);
Ts = 0.005;
U_max = 5;
T_max = 0.04;
normalization_impulse = 1 / (U_max * T_max);
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

% Frekvencni charakteristika
data_shaft_freq = readmatrix("./data/hridel_frek_1.csv"); 
data_shaft_freq_cleaned = data_shaft_freq(:,[4:7]);
%%
Ts = 0.005;
%x = linspace(0, length(data_shaft_freq_cleaned(2000:3000)) * Ts, length(data_shaft_freq_cleaned(2000:3000)));

figure 
hold on
plot((2057:2257)*Ts-2057*Ts, data_shaft_freq_cleaned(2057:2257, 3)) % Budici napeti
plot((2057:2257)*Ts-2057*Ts, data_shaft_freq_cleaned(2057:2257, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika (Amp = 1; f = 1 Hz)")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure 
hold on
plot((6400:6458)*Ts-6400*Ts, data_shaft_freq_cleaned(6400:6458, 3)) % Budici napeti
plot((6400:6458)*Ts-6400*Ts, data_shaft_freq_cleaned(6400:6458, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika (Amp = 2; f = 3.5 Hz)")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure 
hold on
plot((14269:14320)*Ts-14269*Ts, data_shaft_freq_cleaned(14269:14320, 3)) % Budici napeti
plot((14269:14320)*Ts-14269*Ts, data_shaft_freq_cleaned(14269:14320, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika (Amp = 1.5; f = 4 Hz)")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Rychlost hridele setrvacniku")

%%
% SC2FA
numerator = [-8.1177, 332.74];
denominator = [1, 13.025, 791.44];

sys = tf(numerator, denominator);

data_shaft_sc2fa = readmatrix("./data/hridel_sc2fa_1.csv");
data_frek = data_shaft_sc2fa(:,5:6);
data_nyq = data_frek(find(data_shaft_sc2fa(:,7)==1),:);

figure
hold on
plot(data_nyq(:,1), data_nyq(:,2))
%nyquist(sys)
title("Nyquist diagram")
xlabel("Re")
ylabel("Im")
%grid on


% Porovnani prechodovych charakteristik
%Ts=0.01;
[s,t]=step(sys,0:Ts:(length(data_shaft_step_cleaned_normalized(1940:2200, 3))-1)*Ts);
figure
hold on
plot(0:Ts:(length(data_shaft_step_cleaned_normalized(1940:2200, 3))-1)*Ts, data_shaft_step_cleaned_normalized(1940:2200, 2))
plot(t,s)
title("Porovnani prechodovych charakteristik")
xlim([0 2])
xlabel("t [s]")
ylabel("y (t)")


%% Pruzny pas
% Prechodova charakteristika
data_laser_belt_step = readmatrix("./data/pruzny_pas_step_b1.csv");
data_induction_belt_step = readmatrix("./data/pruzny_pas_step_b2.csv"); 

data_belt_step_cleaned = data_induction_belt_step(:,[4:7]);

Ts = 0.01;
x = linspace(0, length(data_belt_step_cleaned) * Ts, length(data_belt_step_cleaned));
U_max = 10;
normalization_step = 1 / U_max;

data_belt_step_cleaned_normalized = data_belt_step_cleaned * normalization_step;

figure
hold on
plot(x, data_belt_step_cleaned_normalized(:, 4)) 
plot(x, data_belt_step_cleaned_normalized(:, 1))
%plot(x, data_belt_step_cleaned_normalized(:, 3))
plot(x, data_belt_step_cleaned_normalized(:, 2))
%xlim([0 100])
title("Prechodova charakteristika")
xlabel("t [s]")
%legend("Budici napeti", "Rychlost hridele motoru", "Rychlost hridele setrvacniku")


% Frekvencni charakteristika
data_belt_freq = readmatrix("./data/pruzny_pas_freq_c2.csv");
data_belt_freq_cleaned = data_belt_freq(:,[4:7]);

x = linspace(0, length(data_belt_freq_cleaned) * Ts, length(data_belt_freq_cleaned));

figure
hold on
plot(x, data_belt_freq_cleaned(:, 4)) 
plot(x, data_belt_freq_cleaned(:, 1))
%plot(x, data_belt_freq_cleaned(:, 3))
plot(x, data_belt_freq_cleaned(:, 2))
%xlim([0 100])
title("Frekvencni charakteristika")
xlabel("t [s]")
xlim([0 100])


% SC2FA
numerator = [-40.198, 851.79];
denominator = [1, 52.736, 168.38];

sys = tf(numerator, denominator);

data_belt_sc2fa = readmatrix("./data/pruzny_pas_sc2fa_1.csv");
data_frek = data_belt_sc2fa(:,5:6);
data_nyq = data_frek(find(data_belt_sc2fa(:,7)==1),:);

x = linspace(0, length(data_nyq) * Ts, length(data_nyq)); 
x_frek = linspace(0, length(data_belt_sc2fa)*Ts, length(data_belt_sc2fa)); 

figure
hold on
plot(x, data_nyq(:,1))
plot(x, data_nyq(:,2))
plot(x_frek, data_belt_sc2fa(:,4))
%nyquist(sys)
%title("Nyquist diagram")
xlabel("t [s]")
%ylabel("Im")
%grid on

%% Teplomer
data_thermometer = readmatrix("./data/teplomer_data_all.csv"); 
data_thermometer_cleaned = data_thermometer(:,[3:5]); % without NaN values [..., temperature, ...]
data_therm2=data_thermometer_cleaned(:,3);
Ts = 0.1; % perioda vzorkovani

% Zavislost teploty na napeti
[max_temp, max_idx] = max(data_thermometer_cleaned(:, 2)); 
[min_temp, min_idx] = min(data_thermometer_cleaned(:, 2));

%normalization = 1 / (max_temp - min_temp);

data_thermometer_90_to_25_temperature = data_thermometer_cleaned(3460:max_idx+500, 2)/70 -data_thermometer_cleaned(3460, 2)/70;
data_thermometer_90_to_25_voltage = data_thermometer_cleaned(3460:max_idx+500, 3)/10;

x = linspace(0, length(data_thermometer_90_to_25_temperature) * Ts, length(data_thermometer_90_to_25_temperature));
figure
hold on
plot(x, data_thermometer_90_to_25_temperature)
plot(x, data_thermometer_90_to_25_voltage)
xlabel("t [s]")
ylabel("y_1(t), y_2(t)")
title("Prechodova charakteristika")

%% Eddy current
% Prechodova charakteristika
data_eddy_step = readmatrix("./data/Eddy_step.csv"); 
data_eddy_step_cleaned = data_eddy_step(:,[4:7]);
Ts = 0.001;
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
xlabel("t [s]")
ylabel("u(t), y(t)")

% Frevencni charakteristika
data_eddy_freq = readmatrix("./data/Eddy_frek_14112023_new.csv"); 
data_eddy_freq_cleaned = data_eddy_freq(:,[4:7]);

x = linspace(0, length(data_eddy_freq_cleaned) * Ts, length(data_eddy_freq_cleaned));

figure
hold on
plot(x, data_eddy_freq_cleaned(:, 3))
plot(x, data_eddy_freq_cleaned(:, 2))
plot(x, data_eddy_freq_cleaned(:, 1))
title("Frekvencni charakteristika")
xlabel("t")

% Blok FRID
data_eddy_FRID = readmatrix("./data/Eddy_FRID_freq.csv");
data_eddy_FRID_cleaned = data_eddy_FRID(1:23924,[4:7]);

x = linspace(0, length(data_eddy_FRID_cleaned) * Ts, length(data_eddy_FRID_cleaned));

figure
hold on
plot(data_eddy_FRID_cleaned(:, 2), data_eddy_FRID_cleaned(:, 3))
title("Nyquist diagram")
xlabel("Re")
ylabel("Im")


figure
plot(x, data_eddy_FRID_cleaned(:, 1))
title("Eddy current (FRID)")
ylabel("Frekvence")