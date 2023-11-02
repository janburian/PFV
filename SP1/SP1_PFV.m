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
y_manufacturer_ultra = (4/350) * x - 1.42

figure;
scatter(distances, measured_values_1, '*', "blue")
hold on
scatter(distances, measured_values_2, '+', "red")
hold on
plot(x, y); 
hold on
plot(x, y_manufacturer_ultra)
xlabel("d [mm]")
ylabel("U [V]")
title("Staticka charakteristika")
legend("Mereni 1", "Mereni 2", "Aproximacni primka", "Staticka charakteristika dana vyrobcem")

% Chyba opakovatelnosti
data_ultra_repetition = readmatrix("./data/Ultra_opak60cm.csv");
data_ultra_repetition_voltage = data_ultra_repetition(:, 4);
%data_ultra_repetition_voltage_avg = sum(data_ultra_repetition(:, 4)) / length(data_ultra_repetition_voltage);
p_2 = polyfit([measured_values_1, measured_values_2], [distances, distances], 1);
y = p_2(1) * data_ultra_repetition_voltage + p_2(2);
 
[muV,varV,deltaV,dV]=opak(y,300,1000)

%% Mereni teploty
data_thermometer = readmatrix("./data/teplomer_data_all.csv"); 
data_thermometer_cleaned = data_thermometer(:,[3:5]); % without NaN values [..., temperature, ...]
data_therm2=data_thermometer_cleaned(:,3);
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

data_thermometer_90_to_25_temperature = [data_thermometer_cleaned(max_idx:end, 2); flip(data_thermometer_cleaned(1:500, 2))];
data_thermometer_90_to_25_voltage = [data_thermometer_cleaned(max_idx:end, 3); flip(data_thermometer_cleaned(1:500, 3))];

% p_3 = polyfit(data_thermometer_cleaned(max_idx:end, 2), data_thermometer_cleaned(max_idx:end, 3), 3);
% p_4 = polyfit(data_thermometer_cleaned(max_idx:end, 3), data_thermometer_cleaned(max_idx:end, 2), 3);

p_3 = polyfit(data_thermometer_90_to_25_temperature, data_thermometer_90_to_25_voltage, 3);

x = linspace(25, 94, 1000); % Adapt n for resolution of graph
y = p_3(1) * x.^3 + p_3(2) * x.^2 + p_3(3) * x + p_3(4);


figure
plot(data_thermometer_90_to_25_temperature, data_thermometer_90_to_25_voltage)
hold on
plot(x, y)
xlabel("Teplota namerena referencnim snimacem [°C]")
ylabel("Napeti namerene polovodicovym snimacem U [V]")
title("Zavislost teploty na napeti")
legend("Staticka charakteristika", "Aproximacni polynom")

y_inv=p_4(1) * data_therm2.^3 + p_4(2) * data_therm2.^2 + p_4(3) * data_therm2 + p_4(4);
figure
plot(linspace(0, length(data_thermometer_cleaned) * Ts, length(data_thermometer_cleaned)),y_inv)
xlabel("Cas [s]")
ylabel("Teplota [°C]")
title("Vystup referencniho a polovodicoveho snimace")
yline(60, '--', 'Odkryti vodni lazne');
legend("Teplota namerena polovodicovym snimacem T [°C]")

% Opak chyba
y_inv=y_inv(1:2000);
[muT,varT,deltaT,dT]=opak(y_inv,0,100)

%% Elektromechanicka soustava modelu pruzne hridele
data_shaft = readmatrix("./data/hridel_ukol_c.csv"); 
data_shaft_cleaned = data_shaft(:,[4:6]); % without NaN and 0 values
Ts = 0.02; % perioda vzorkovani

x = linspace(0, length(data_shaft_cleaned) * Ts, length(data_shaft_cleaned));
figure
plot(x, data_shaft_cleaned(:, 1), "-")
hold on 
plot(x, data_shaft_cleaned(:, 2), "-")
hold on 
plot(x, data_shaft_cleaned(:, 3), "-")
%xlim([180 400])
xlabel("Cas [s]")
ylabel("U [V]")
title("Casovy prubeh signalu")
legend("IRC1", "IRC2", "Vstupni signal")

% Staticka charakteristika
voltage_range = find(data_shaft_cleaned(:,3) > -9 & data_shaft_cleaned(:,3) < 9);
p_6 = polyfit(data_shaft_cleaned(voltage_range,3), data_shaft_cleaned(voltage_range,1), 1);
x = linspace(-9, 9, 1000); % Adapt n for resolution of graph
y = p_6(1) * x + p_6(2);

figure
plot(data_shaft_cleaned(:, 3), data_shaft_cleaned(:, 1));
hold on
plot(x, y, 'LineWidth', 1.3)
xlabel("Napeti na motor [V]")
ylabel("Uhlova rychlost motoru [rad\cdots^{-1}]")
title("Staticka charakteristika")
legend("Staticka charakteristika", "Aproximacni primka")

% 
data_shaft_rotation = readmatrix("./data/hridel_10_otaceni.csv"); 
data_shaft_rotation_cleaned = data_shaft(:,[4:5]); % without NaN and 0 values
Ts = 0.02; % perioda vzorkovani

%% Pruzny pas
data_belt = readmatrix("./data/pruzny_pas_C-a.csv"); 
data_belt_cleaned = data_belt(:,[6:7]); % without NaN and 0 values

voltage_range = find(data_belt_cleaned(:,2) > -7 & data_belt_cleaned(:,2) < 7);
p_5 = polyfit(data_belt_cleaned(voltage_range,1), data_belt_cleaned(voltage_range,2), 1);
p_6 = polyfit(data_belt_cleaned(voltage_range,2), data_belt_cleaned(voltage_range,1), 1);

x = linspace(90, 115, 1000); % Adapt n for resolution of graph
y = p_5(1) * x + p_5(2);

Ts = 0.01;

figure
plot(data_belt_cleaned(:,1), data_belt_cleaned(:,2))
hold on
plot(x, y)
title("Staticka charakteristika")
xlabel("Laserovy snimac d [mm]")
ylabel("Indukcni snimac U [V]")
legend("Staticka charakteristika", "Primkova aproximace")

figure
x = linspace(0, length(data_belt_cleaned) * Ts, length(data_belt_cleaned));
plot(x, data_belt_cleaned(:,1))
title("Vzdalenost volne kladky od laseroveho snimace")
xlabel("Cas t [s]")
ylabel("Vzdalenost d [mm]")

% figure
% x = linspace(0, length(data_belt_cleaned) * Ts, length(data_belt_cleaned));
% plot(data_belt_cleaned(:,2), data_belt_cleaned(:,1))
%%

%data_belt_b = readmatrix("./data/pruzny_pas_b_Data_correct.csv"); 
data_belt_b = readmatrix("./data/pruzny_pas_b_Data_correct.csv"); 
p_6b=polyfit(data_belt_b(:,7),(p_6(1)*data_belt_b(:,4)+p_6(2)),1);

LinChar=p_6b(1).*(-10:10) +p_6b(2);
figure
hold on
plot( data_belt_b(:,7),(p_6(1)*data_belt_b(:,4)+p_6(2)))
plot(-10:10, LinChar)
title("Staticka charakteristika soustavy")
xlabel("U [V]")
ylabel("d [mm]")
legend("Staticka charakteristika", "Aproximacni primka");
%%
data_belt_c = readmatrix("./data/pruzny_pas_kontrola.csv"); 
data_belt_c_cleaned = data_belt_c(1:15538, [4:7]);%155538

x=-1.5:0.1:1.5;
y=polyfit(data_belt_c_cleaned(:,4),(3.02*data_belt_c_cleaned(:,3)+101.972),1);
y=y(1)*x+y(2);
figure
hold on
plot( data_belt_c_cleaned(:,4),(3.02*data_belt_c_cleaned(:,3)+101.972));
plot(x,y)
plot([0,0], [101.5,110], 'k--')
title("Staticka charakteristika soustavy")
xlabel("U [V]")
ylabel("d [mm]")
legend("Staticka charakteristika", "Aproximacni primka");


voltage_range = find(data_belt_c_cleaned(:,4) >-0.7 & data_belt_c_cleaned(:,4) < 1.18);
LinChar=polyfit(data_belt_c_cleaned(voltage_range,4),data_belt_c_cleaned(voltage_range,1) - data_belt_c_cleaned(voltage_range,2),1);

x=-0.7:0.1:1.2;
y=LinChar(1)*x+LinChar(2);
figure
hold on


plot( data_belt_c_cleaned(:,4),data_belt_c_cleaned(:,1) - data_belt_c_cleaned(:,2))
plot(x,y, "LineWidth",1.3)
title("Rozdil rychlosti motoru v zavislosti na budicim signalu")
xlabel("Rozdil v budicim signalu [V]")
ylabel("Rozdil rychlosti motoru [m\cdot s^{-1}]")
legend("Staticka charakteristika prokluzu", "Aproximacni primka");

%%
data_belt_chyba = readmatrix("./data/pruzny_pas_chyba_opak.csv");
chyba_laser= data_belt_chyba(:,6);
chyba_ind= p_6(1)*data_belt_chyba(:,7)+p_6(2);
[mLaser, varLaser, DLaser, dLaser]=opak(chyba_laser,88.2, 125.6)
[mInd, varInd, DInd, dInd]=opak(chyba_ind,88.2, 125.6)

%% Indukcni snimac vzdalenosti typu PR6423 (Eddy current)
distances_eddy = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]; % [mm]
coeff = 3.3 / 13.3; 
%measured_values_eddy = -1/coeff * [3.44, 3.14, 2.88, 2.62, 2.35, 2.12, 1.9, 1.75, 1.55, 1.29, 1.03]; % [V]
measured_values_eddy = [3.44, 3.14, 2.88, 2.62, 2.35, 2.12, 1.9, 1.75, 1.55, 1.29, 1.03]; % [V]

p_eddy = polyfit(distances_eddy, measured_values_eddy, 1);
p_eddy_2 = polyfit( measured_values_eddy,distances_eddy, 1);
x = linspace(0.5, 2.5, 11); % Adapt n for resolution of graph
y_eddy = p_eddy(1) * x + p_eddy(2);
y_eddy_manufacturer = -coeff*linspace(-18, -2, 11);

p_eddy_prokluz=polyfit( x,y_eddy_manufacturer, 1);

figure;
scatter(distances_eddy,measured_values_eddy, '*', "blue")
hold on
plot(x, y_eddy)
hold on 
plot(x, y_eddy_manufacturer)
xlabel("d [mm]")
ylabel("U [V]")
title("Staticka charakteristika indukcniho snimace vzdalenosti")
legend("Namerena data", "Aproximacni primka", "Staticka charakteristika dana vyrobcem")
%%
data_eddy = readmatrix("./data/eddy_5.csv");
data_eddy_opak=p_eddy_2(1)*data_eddy(:,4)+p_eddy_2(2);
[mEddy,varEddy, DeltaEddy, dEddy]=opak(data_eddy_opak,0.5,2.5)

%% Chyby opakovatelnosti
function [me,vr,Deltax, dx] = opak(data,dmin,dmax)
me=mean(data);
vr=var(data);
Deltax=2*sqrt(vr);
dx=Deltax/(dmax-dmin)*100;
fprintf("Charakteristiky chyb opakovatelnosti:")
end
