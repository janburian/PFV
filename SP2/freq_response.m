clc
clear all
close all

%% Frequency response
%% Shaft
% 1
max_y_1 = 0.6;
min_y_1 = -0.61;
max_u_1 = 1;
min_u_1 = -1;
delta_t_1 = 0.06;
f_1 = 1;

Amp_1 = count_amplitude(max_y_1, min_y_1, max_u_1, min_u_1);
Amp_dB_1 = count_amplitude_dB(Amp_1)
phi_1 = count_phase(f_1, delta_t_1)


% 2
max_y_2 = 1.682; 
min_y_2 = -1.722; 
max_u_2 = 2;
min_u_2 = -2;
delta_t_2 = 0.0550;
f_2 = 3.5; 

Amp_2 = count_amplitude(max_y_2, min_y_2, max_u_2, min_u_2);
Amp_dB_2 = count_amplitude_dB(Amp_2)
phi_2 = count_phase(f_2, delta_t_2)


% 3
max_y_3 = 1.481;
min_y_3 = -1.451;
max_u_3 = 1.5;
min_u_3 = -1.5;
delta_t_3 = 0.0550;
f_3 = 4;

Amp_3 = count_amplitude(max_y_3, min_y_3, max_u_3, min_u_3);
Amp_dB_3 = count_amplitude_dB(Amp_3)
phi_3 = count_phase(f_3, delta_t_3)

%% Belt
% 1
max_y_1 = 0;
min_y_1 = -2.443;
max_u_1 = 1;
min_u_1 = -1;
delta_t_1 = 0.24;
f_1 = 1;

Amp_1 = count_amplitude(max_y_1, min_y_1, max_u_1, min_u_1);
Amp_dB_1 = count_amplitude_dB(Amp_1)
phi_1 = count_phase(f_1, delta_t_1)


% 2
max_y_2 = 0.121; 
min_y_2 = -2.204; 
max_u_2 = 2;
min_u_2 = -2;
delta_t_2 = 0.14;
f_2 = 3.5; 

Amp_2 = count_amplitude(max_y_2, min_y_2, max_u_2, min_u_2);
Amp_dB_2 = count_amplitude_dB(Amp_2)
phi_2 = count_phase(f_2, delta_t_2)


% 3
max_y_3 = 0.198;
min_y_3 = -2.316;
max_u_3 = 3;
min_u_3 = -3;
delta_t_3 = 0.1;
f_3 = 6;

Amp_3 = count_amplitude(max_y_3, min_y_3, max_u_3, min_u_3);
Amp_dB_3 = count_amplitude_dB(Amp_3)
phi_3 = count_phase(f_3, delta_t_3)

%% Eddy current
% 1
max_y_1 = 0.1465;
min_y_1 = -0.1026;
max_u_1 = 0.1;
min_u_1 = -0.1;
delta_t_1 = 0.052;
f_1 = 10;

Amp_1 = count_amplitude(max_y_1, min_y_1, max_u_1, min_u_1);
Amp_dB_1 = count_amplitude_dB(Amp_1)
phi_1 = count_phase(f_1, delta_t_1)


% 2
max_y_2 = 0.8449; 
min_y_2 = -0.7937; 
max_u_2 = 0.2;
min_u_2 = -0.2;
delta_t_2 = 0.016;
f_2 = 50; 

Amp_2 = count_amplitude(max_y_2, min_y_2, max_u_2, min_u_2);
Amp_dB_2 = count_amplitude_dB(Amp_2)
phi_2 = count_phase(f_2, delta_t_2)

% 3
max_y_3 = 0.193;
min_y_3 = -0.1612;
max_u_3 = 0.3;
min_u_3 = -0.3;
delta_t_3 = 0.03;
f_3 = 100;

Amp_3 = count_amplitude(max_y_3, min_y_3, max_u_3, min_u_3);
Amp_dB_3 = count_amplitude_dB(Amp_3)
phi_3 = count_phase(f_3, delta_t_3)

%% Functions
function amplitude = count_amplitude(max_y, min_y, max_u, min_u)
    amplitude = (max_y - min_y) / (max_u - min_u); 
end

function amplitude_dB = count_amplitude_dB(amplitude)
    amplitude_dB = 20 * log10(amplitude); 
end

function phi = count_phase(frequency, delta_t)
    phi = (2 * pi * frequency * delta_t) * (180 / pi); 
    phi = wrapTo180(phi);
end


