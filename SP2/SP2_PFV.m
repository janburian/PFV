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
title("Prechodova charakteristika (hridel)")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Rychlost hridele setrvacniku")

figure
hold on
plot(x, data_shaft_step_cleaned_normalized(:, 3))
plot(x, data_shaft_step_cleaned_normalized(:, 1))
xlim([0 50])
title("Prechodova charakteristika (motor)")
xlabel("t [s]")
ylabel("y(t)")
legend("Budici napeti", "Rychlost hridele motoru")
%{
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

%}
%%
% Sample data points
x_data = x(1:12/Ts);
y_data = data_shaft_step_cleaned_normalized(:, 2);
y_data=y_data(1:12/Ts);
% Point at which to draw the tangent
%ipt = findchangepts(data_shaft_step_cleaned_normalized(1:2100, 2));
given_point_index = findchangepts(y_data); % Adjust as needed

% Define a neighborhood around the given point
neighborhood_size = 2;
neighborhood_indices = given_point_index - neighborhood_size:given_point_index + neighborhood_size;

% Fit a linear function to the neighborhood
coefficients = polyfit(x_data(neighborhood_indices), y_data(neighborhood_indices)', 1);
tangent_line = polyval(coefficients, x_data);

% Plot the data and the tangent line
figure;
hold on;
plot(x,data_shaft_step_cleaned_normalized(:, 3))
plot(x_data, y_data, 'DisplayName', 'Data');
[v,i]=max(y_data);
scatter((i-1)*Ts,v,'gx');
scatter(9.709,0,'bo');
scatter(10.05,0.694,'ko');
yline(mean(y_data(end-100:end))*0.95,'g--')
yline(mean(y_data(end-100:end))*1.05,'g--')
plot(x_data, tangent_line,'m');
xlim([9.5 10.5])
ylim([-0.2 1.1])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Vystup senzoru","\sigma_{max}","D","T_u","\sigma_{5%}")
grid on;
hold off;



%%
% Impulsni charakteristika
data_shaft_impulse = readmatrix("./data/hridel_imp_1.csv"); 
data_shaft_impulse_cleaned = data_shaft_impulse(:,[4:7]);
Ts = 0.005;
U_max = 5;
T_max = 0.03;
normalization_impulse = 1 / (U_max * T_max);
x = linspace(0, length(data_shaft_impulse_cleaned) * Ts, length(data_shaft_impulse_cleaned));

data_shaft_impulse_cleaned_normalized = data_shaft_impulse_cleaned * normalization_impulse;

figure
hold on
plot(x, data_shaft_impulse_cleaned_normalized(:, 3))
%plot(x, data_shaft_impulse_cleaned_normalized(:, 1))
plot(x, data_shaft_impulse_cleaned_normalized(:, 2))
title("Impulsni charakteristika")
ylim([-0.9 1.7])
xlim([10 30])
xlabel("t[s]")
ylabel("u(t),y(t)")
legend("Budici napeti", "Rychlost hridele setrvacniku")


%%
% Sample data points
x_data = x(1:13.5/Ts);
y_data = data_shaft_impulse_cleaned_normalized(:, 2);
y_data=y_data(1:13.5/Ts);
% Point at which to draw the tangent
%ipt = findchangepts(data_shaft_step_cleaned_normalized(1:2100, 2));
given_point_index = findchangepts(y_data); % Adjust as needed

% Define a neighborhood around the given point
neighborhood_size = 2;
neighborhood_indices = given_point_index - neighborhood_size:given_point_index + neighborhood_size;

% Fit a linear function to the neighborhood
coefficients = polyfit(x_data(neighborhood_indices), y_data(neighborhood_indices)', 1);
tangent_line = polyval(coefficients, x_data);

% Plot the data and the tangent line
figure;
hold on;
plot(x,data_shaft_impulse_cleaned_normalized(:, 3))
plot(x_data, y_data, 'DisplayName', 'Data');
[v,i]=max(y_data);
scatter((i-1)*Ts,v,'gx');
scatter(12.17,0,'bo');
scatter(12.55,1.049,'ko');
yline(mean(y_data(end-100:end))*0.95,'g--')
yline(mean(y_data(end-100:end))*1.05,'g--')
plot(x_data, tangent_line,'m');
xlim([11.5 13.5])
ylim([-0.2 1.8])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Vystup senzoru","\sigma_{max}","D","T_u","\sigma_{5%}")
grid on;
hold off;







%% 
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

l=59944;
figure
hold on
plot(data_nyq(401:l,1), data_nyq(401:l,2))
plot(data_nyq(l+1:end,1), data_nyq(l+1:end,2))
%nyquist(sys)
title("Nyquistuv diagram")
xlabel("Re")
ylabel("Im")
legend("Namerena data", "Krivka odhadnuta blokem sc2fa")
%grid on

%% Bode diagram
data_freq = data_shaft_sc2fa(:,4)*10;
omega = data_freq(find(data_shaft_sc2fa(:,7)==1),:);
magnitudes_measurements = sqrt(data_nyq(401:l,1).^2 + data_nyq(401:l,2) .^2);
phases_measurements = atan2(data_nyq(401:l,2), data_nyq(401:l,1)) * (180 / pi); 
omega = omega(401:l,1);
%omega = logspace(-2, 3, length(magnitudes_measurements)); % Frequency vector in logarithmic scale TODO: change this
% Bode plot
figure;

% Magnitude plot
subplot(2, 1, 1);
hold on
semilogx(omega, 20 * log10(magnitudes_measurements),'k');
bodemag(sys,omega)
title('Bode Plot - Magnitude');
xlabel('Frequency ');
ylabel('Magnitude ');
grid on;

% Phase plot
subplot(2, 1, 2);
hold on

h = bodeplot(sys,omega);

setoptions(h,'MagVisible','off');
semilogx(omega, unwrap(phases_measurements)+360,'k');
ylim([-180 360])
title('Bode Plot - Phase');
xlabel('Frequency ');
ylabel('Phase ');
grid on;
%%

% Porovnani prechodovych charakteristik
%Ts=0.01;
[s,t]=step(sys,0:Ts:(length(data_shaft_step_cleaned_normalized(1940:2200, 3))-1)*Ts);
figure
hold on
plot(0:Ts:(length(data_shaft_step_cleaned_normalized(1940:2200, 3))-1)*Ts, data_shaft_step_cleaned_normalized(1940:2200, 2))
plot(t,s)
title("Porovnani prechodovych charakteristik")
legend("Namerena data", "Identifikovana data")
xlabel("t [s]")
ylabel("y (t)")


%% Pruzny pas
% Prechodova charakteristika
data_laser_belt_step = readmatrix("./data/pruzny_pas_step_b1.csv");
data_induction_belt_step = readmatrix("./data/pruzny_pas_step_b2.csv"); 

data_belt_step_cleaned = data_laser_belt_step(:,[4:7]);

Ts = 0.01;
x = linspace(0, length(data_belt_step_cleaned) * Ts, length(data_belt_step_cleaned));
U_max = 1;
normalization_step = 1 / U_max;

data_belt_step_cleaned_normalized = data_belt_step_cleaned * normalization_step;

figure
hold on
plot(x, data_belt_step_cleaned_normalized(:, 4)) 
%plot(x, data_belt_step_cleaned_normalized(:, 1))
plot(x, data_belt_step_cleaned_normalized(:, 3)-109.6)
%plot(x, data_belt_step_cleaned_normalized(:, 2))
%xlim([0 100])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Vzdalenost kladky")
%legend("Budici napeti", "Rychlost hridele motoru", "Rychlost hridele setrvacniku")

%%
% Sample data points
x_data = x(1:17/Ts);
y_data = data_belt_step_cleaned_normalized(:, 3)-109.6;
y_data=y_data(1:17/Ts);
% Point at which to draw the tangent
%ipt = findchangepts(data_shaft_step_cleaned_normalized(1:2100, 2));
given_point_index = findchangepts(y_data)-3; % Adjust as needed

% Define a neighborhood around the given point
neighborhood_size = 2;
neighborhood_indices = given_point_index - neighborhood_size:given_point_index + neighborhood_size;

% Fit a linear function to the neighborhood
coefficients = polyfit(x_data(neighborhood_indices), y_data(neighborhood_indices)', 1);
tangent_line = polyval(coefficients, x_data);

% Plot the data and the tangent line
figure;
hold on;
plot(x, data_belt_step_cleaned_normalized(:, 4))
plot(x_data, y_data, 'DisplayName', 'Data');
[v,i]=max(y_data);
scatter((i-1)*Ts,v,'gx');
scatter(13.885,0,'bo');
scatter(15.76,4.198,'ko');
yline(mean(y_data(end-100:end))*0.95,'g--')
yline(mean(y_data(end-100:end))*1.05,'g--')
plot(x_data, tangent_line,'m');
xlim([13 17])
ylim([-0.2 5])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Vystup senzoru","\sigma_{max}","D","T_u","\sigma_{5%}")
grid on;
hold off;



%% Frekvencni charakteristika
data_belt_freq = readmatrix("./data/pruzny_pas_freq_c2.csv");
data_belt_freq_cleaned = data_belt_freq(:,[4:7]);

x = linspace(0, length(data_belt_freq_cleaned) * Ts, length(data_belt_freq_cleaned));

figure
hold on
%plot(x, data_belt_freq_cleaned(:, 4)) 
plot(x, data_belt_freq_cleaned(:, 1)-110.36)
%plot(x, data_belt_freq_cleaned(:, 3))
plot(x, data_belt_freq_cleaned(:, 2))
%xlim([0 100])
title("Frekvencni charakteristika")
legend("Budici napeti", "Odezva")
ylabel("u(t), y(t)")
xlabel("t [s]")
xlim([0 100])


% SC2FA
numerator = [-40.198, 851.79];
denominator = [1, 52.736, 168.38];

sys = tf(numerator, denominator);

data_belt_sc2fa = readmatrix("./data/pruzny_pas_sc2fa_1.csv");
data_frek = data_belt_sc2fa(:,4:6);
data_nyq = data_frek(find(data_belt_sc2fa(:,7)==1),2:3);

frek_bode=data_frek(find(data_belt_sc2fa(:,7)==1),1);
%x = linspace(0, length(data_nyq) * Ts, length(data_nyq)); 
x_frek = linspace(0, length(data_belt_sc2fa)*Ts, length(data_belt_sc2fa)); 

figure
hold on
plot(0:Ts:100249*Ts, data_nyq(1:100250,1))
plot(0:Ts:100249*Ts, data_nyq(1:100250,2))
plot(0:Ts:100249*Ts, data_belt_sc2fa(1:100250,4))
%nyquist(sys)
title("Mereni frekvencni charakteristiky pomoci sc2fa bloku")
legend("Re", "Im", "Frekvence")
xlabel("t [s]")

%ylabel("Im")
%grid on
%%
figure
hold on
plot(data_nyq(601:100250,1),data_nyq(601:100250,2))
plot(data_nyq(102100:end,1),data_nyq(102100:end,2))
legend("Namerena data", "Krivka odhadnuta blokem sc2fa")
xlabel("Re")
ylabel("Im")

%% Bode diagram
data_freq = frek_bode(601:100250);
omega = data_freq;
magnitudes_measurements = sqrt(data_nyq(601:100250,2).^2 + data_nyq(601:100250,1) .^2);
phases_measurements = atan2(data_nyq(601:100250,2), data_nyq(601:100250,1)) * (180 / pi); 
%omega = omega(401:l,1);
%omega = logspace(-2, 3, length(magnitudes_measurements)); % Frequency vector in logarithmic scale TODO: change this
% Bode plot
figure;

% Magnitude plot
subplot(2, 1, 1);
hold on
semilogx(omega, 20 * log10(magnitudes_measurements),'k');
bodemag(sys,omega)
scatter(1,1.7379,'xm');
scatter(3.5,-4.7127,'xm');
scatter(6,-7.5557,'xm');
title('Bode Plot - Magnitude');
xlabel('Frequency ');
ylabel('Magnitude ');
grid on;

% Phase plot
subplot(2, 1, 2);
hold on

h = bodeplot(sys,omega);

setoptions(h,'MagVisible','off');
semilogx(omega, unwrap(phases_measurements)+360,'k');
scatter(6,216,'xm');
ylim([220 380])
title('Bode Plot - Phase');
xlabel('Frequency ');
ylabel('Phase ');
grid on;
%% Teplomer
data_thermometer = readmatrix("./data/teplomer_data_all.csv"); 
data_thermometer_cleaned = data_thermometer(:,[3:5]); % without NaN values [..., temperature, ...]
data_therm2=data_thermometer_cleaned(:,3);
Ts = 0.1; % perioda vzorkovani

% Zavislost teploty na napeti
[max_temp, max_idx] = max(data_thermometer_cleaned(:, 2)); 
[min_temp, min_idx] = min(data_thermometer_cleaned(:, 2));

%normalization = 1 / (max_temp - min_temp);

data_thermometer_90_to_25_temperature = data_thermometer_cleaned(3460:max_idx+100, 2)/70 -data_thermometer_cleaned(3460, 2)/70;
data_thermometer_90_to_25_voltage = data_thermometer_cleaned(3460:max_idx+100, 3)/10;

x = linspace(0, length(data_thermometer_90_to_25_temperature) * Ts, length(data_thermometer_90_to_25_temperature));
figure
hold on
plot(x, data_thermometer_90_to_25_temperature)
plot(x, data_thermometer_90_to_25_voltage)
xlabel("t [s]")
ylabel("y_1(t), y_2(t)")
title("Prechodova charakteristika")


%%
% Sample data points
x_data = x;
y_data = data_thermometer_90_to_25_temperature;
%y_data=y_data(1:17/Ts);
% Point at which to draw the tangent
%ipt = findchangepts(data_shaft_step_cleaned_normalized(1:2100, 2));
given_point_index = findchangepts(y_data)-45; % Adjust as needed

% Define a neighborhood around the given point
neighborhood_size = 2;
neighborhood_indices = given_point_index - neighborhood_size:given_point_index + neighborhood_size;

% Fit a linear function to the neighborhood
coefficients = polyfit(x_data(neighborhood_indices), y_data(neighborhood_indices)', 1);
tangent_line = polyval(coefficients, x_data);

% Plot the data and the tangent line
figure;
hold on;
plot(x, data_thermometer_90_to_25_voltage)
plot(x_data, y_data, 'DisplayName', 'Data');
[v,i]=max(y_data);
scatter((i-1)*Ts,v,'gx');
scatter(2.205,0,'bo');
scatter(20.45,0.8979,'ko');
yline(mean(y_data(end-100:end))*0.95,'g--')
xline(1.203,'k--')
yline(mean(y_data(end-100:end))*1.05,'g--')
plot(x_data, tangent_line,'m');
%xlim([13 17])
ylim([-0.1 1.1])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Vystup senzoru","\sigma_{max}","D","T_u","\sigma_{5%}")
grid on;
hold off;






%% Eddy current
% model
num=1.0e+06 *[0.1781 4.4089];
den=[1 91.36 1.507e05 3.158e06];
d=0.003;
model=tf(num,den,'Inputdelay',d);
% Prechodova charakteristika
data_eddy_step = readmatrix("./data/Eddy_step.csv"); 
data_eddy_step_cleaned = data_eddy_step(:,4:7);
Ts = 0.001;
x = linspace(0, length(data_eddy_step_cleaned) * Ts, length(data_eddy_step_cleaned));
U_max = 5;
normalization_step = 1 / U_max;
%{
figure
hold on
plot(x, (data_eddy_step_cleaned(:, 3)-0.5)*5)
%plot(x, data_eddy_step_cleaned(:, 2))
plot(x, -(data_eddy_step_cleaned(:, 1)-2.27)*5)
xlim([0 1.5])
ylim([-2.2 2.2])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Vystup senzoru")
%}

% Sample data points
x_data = x(1:190);
y_data = -(data_eddy_step_cleaned(:, 1)-2.27)*5;
y_data=y_data(1:190);
% Point at which to draw the tangent
%ipt = findchangepts(data_shaft_step_cleaned_normalized(1:2100, 2));
given_point_index = findchangepts(y_data); % Adjust as needed

% Define a neighborhood around the given point
neighborhood_size = 2;
neighborhood_indices = given_point_index - neighborhood_size:given_point_index + neighborhood_size;

% Fit a linear function to the neighborhood
coefficients = polyfit(x_data(neighborhood_indices), y_data(neighborhood_indices)', 1);
tangent_line = polyval(coefficients, x_data);

% Plot the data and the tangent line
figure;
hold on;
plot(x, (data_eddy_step_cleaned(:, 3)-0.5)*5)
plot(x_data, y_data, 'DisplayName', 'Data');
[v,i]=max(y_data);
scatter((i-1)*Ts,v,'gx');
scatter(0.0085,0,'bo');
scatter(0.097,1.435,'ko');
yline(mean(y_data(end-100:end))*0.95,'g--')
yline(mean(y_data(end-100:end))*1.05,'g--')
plot(x_data, tangent_line,'m');
xlim([0 0.13])
ylim([-0.2 2.2])
title("Prechodova charakteristika")
xlabel("t [s]")
ylabel("u(t), y(t)")
legend("Budici napeti", "Vystup senzoru","\sigma_{max}","D","T_u","\sigma_{5%}")
grid on;
hold off;

%% porovnani prechodu
[a,b]=step(model);
figure
hold on
plot(0:Ts:(length(data_eddy_step_cleaned(0.005/Ts:end, 3))-1)*Ts, (data_eddy_step_cleaned(0.005/Ts:end, 3)-0.5)*5, 'm')
plot(x_data(1:end-0.006/Ts+1), y_data(0.006/Ts:end));
%step(model)
plot(b,a,'b--')
xlim([0 0.15])
ylim([-0.2 2.2])
ylabel("u(t),y(t)")
xlabel("t[s]")
title("Porovnani namerene prechodove charakteristiky a charakteristiky modelu")
legend("Vstupni signal", "Namerena data", "Vystup modelu")
%% Frevencni charakteristika
data_eddy_freq = readmatrix("./data/Eddy_frek_14112023_new.csv"); 
data_eddy_freq_cleaned = data_eddy_freq(:,[4:7]);



x = linspace(0, length(data_eddy_freq_cleaned) * Ts, length(data_eddy_freq_cleaned));

figure
hold on
plot(x, data_eddy_freq_cleaned(:, 3)-0.5)
%plot(x, data_eddy_freq_cleaned(:, 2)+5)
plot(x, data_eddy_freq_cleaned(:, 1)-2)
title("Frekvencni charakteristika")
legend("Budici napeti", "Odezva")
ylabel("u(t), y(t)")
xlabel("t [s]")



%% Blok FRID

[re,im,w]=nyquist(model);
[mag,phase,wout] = bode(model);
data_eddy_FRID = readmatrix("./data/Eddy_FRID_freq.csv");
data_eddy_FRID_cleaned = data_eddy_FRID(1:23924,[4:7]);

x = linspace(0, length(data_eddy_FRID_cleaned) * Ts, length(data_eddy_FRID_cleaned));
n=2000;

figure
hold on
plot(data_eddy_FRID_cleaned(n:end, 2), data_eddy_FRID_cleaned(n:end, 3))
title("Nyquistuv diagram")
xlabel("Re")
ylabel("Im")

figure
hold on
plot(data_eddy_FRID_cleaned(n:end, 2), data_eddy_FRID_cleaned(n:end, 3))
plot(squeeze(re),squeeze(im))
title("Nyquistuv diagram")
xlabel("Re")
ylabel("Im")
legend("Namerena data", "Identifikovany model")
%{
figure
plot(x, data_eddy_FRID_cleaned(:, 1))
title("Blok FRID (Eddy current)")
ylabel("Frekvence")
xlabel("t [s]")
%}


%% Bode diagram
data_freq = data_eddy_FRID_cleaned(:,1);
omega = data_freq(n:end);
magnitudes_measurements = sqrt(data_eddy_FRID_cleaned(n:end, 3).^2 + data_eddy_FRID_cleaned(n:end, 2) .^2);
phases_measurements = atan2(data_eddy_FRID_cleaned(n:end, 3), data_eddy_FRID_cleaned(n:end, 2)) * (180 / pi); 
%omega = omega(401:l,1);
%omega = logspace(-2, 3, length(magnitudes_measurements)); % Frequency vector in logarithmic scale TODO: change this
% Bode plot
figure;

% Magnitude plot
subplot(2, 1, 1);
hold on
semilogx(omega, 20 * log10(magnitudes_measurements),'k');
bodemag(sys,omega)
title('Bode Plot - Magnitude');
xlabel('Frequency ');
ylabel('Magnitude ');
grid on;

% Phase plot
subplot(2, 1, 2);
hold on

h = bodeplot(sys,omega);

setoptions(h,'MagVisible','off');
semilogx(omega, unwrap(phases_measurements)+360,'k');
ylim([-180 360])
title('Bode Plot - Phase');
xlabel('Frequency ');
ylabel('Phase ');
grid on;