clc
close all
clear all

%% Hridel 
% Prechodova charakteristika
data_shaft_step = readmatrix("./data/hridel_step_1.csv"); 
data_shaft_step_cleaned = data_shaft_step(:,[4:7]);
Ts = 0.01; % TODO: spravna perioda vzorkovani? DONE: Asi jo!
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
Ts = 0.01; % TODO: spravna perioda vzorkovani? 
U_max = 5;
T_max = 0.04;
normalization_impulse = 1 / (U_max * T_max); % TODO: normalizace
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
%%
Ts=0.005;
x = linspace(0, length(data_shaft_freq_cleaned(2000:3000)) * Ts, length(data_shaft_freq_cleaned(2000:3000)));

figure 
hold on
plot((2057:2257)*Ts-2057*Ts, data_shaft_freq_cleaned(2057:2257, 3)) % Budici napeti
plot((2057:2257)*Ts-2057*Ts, data_shaft_freq_cleaned(2057:2257, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika 1")
%xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure 
hold on
plot((6400:6458)*Ts-6400*Ts, data_shaft_freq_cleaned(6400:6458, 3)) % Budici napeti
plot((6400:6458)*Ts-6400*Ts, data_shaft_freq_cleaned(6400:6458, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika 2")
%xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure 
hold on
plot((14269:14320)*Ts-14269*Ts, data_shaft_freq_cleaned(14269:14320, 3)) % Budici napeti
plot((14269:14320)*Ts-14269*Ts, data_shaft_freq_cleaned(14269:14320, 2)) % Rychlost hridele setrvacniku
title("Frekvencni charakteristika 3")
%xlabel("t")
legend("Budici napeti", "Rychlost hridele setrvacniku")


%%
% SC2FA
numerator = [-8.1177, 332.74];
denominator = [1, 13.025, 791.44];

sys = tf(numerator, denominator)

data_shaft_sc2fa = readmatrix("./data/hridel_sc2fa_1.csv");
data_frek= data_shaft_sc2fa(:,5:6);
data_nyq=data_frek(find(data_shaft_sc2fa(:,7)==1),:);
plot(data_nyq(:,1), data_nyq(:,2))

Ts=0.01;
[s,t]=step(sys,0:Ts:(length(data_shaft_step_cleaned_normalized(1940:2200, 3))-1)*Ts);
figure
hold on
plot(0:Ts:(length(data_shaft_step_cleaned_normalized(1940:2200, 3))-1)*Ts, data_shaft_step_cleaned_normalized(1940:2200, 2))
plot(t,s)
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
data_therm2=data_thermometer_cleaned(:,3);
Ts = 0.1; % perioda vzorkovani

% Zavislost teploty na napeti
[max_temp, max_idx] = max(data_thermometer_cleaned(:, 2)); 
data_thermometer_90_to_25_temperature = data_thermometer_cleaned(3460:max_idx+500, 2)/70 -data_thermometer_cleaned(3460, 2)/70;
data_thermometer_90_to_25_voltage = data_thermometer_cleaned(3460:max_idx+500, 3)/10;


x = linspace(0, length(data_thermometer_90_to_25_temperature) * Ts, length(data_thermometer_90_to_25_temperature));
figure
hold on
plot(x, data_thermometer_90_to_25_temperature)
plot(x, data_thermometer_90_to_25_voltage)
xlabel("t [s]")
ylabel("Napeti namerene polovodicovym snimacem U [V]")
title("Zavislost teploty na napeti")
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

