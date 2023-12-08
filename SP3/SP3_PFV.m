clc 
close all
clear all

%% Shaft - IRC sensor
% frequencies = [2.832, 16.13, 16.77, 1.3719, 8.12, 8.25]; 
% 
% path = "./data/shaft_IRC/"; 
% avg_U_out1 = get_avg_U_out(path, "scope_3.csv");
% avg_U_out2 = get_avg_U_out(path, "scope_4.csv");
% avg_U_out3 = get_avg_U_out(path, "scope_9.csv");
% avg_U_out4 = get_avg_U_out(path, "scope_10.csv");
% %% 
% function [avg_U_out] = get_avg_U_out(path, filename)
%     data = readmatrix(path + filename);
%     data = data(3:end, :);
%     
%     avg_U_out = mean(data(:, 2));
% end

%%
U_out = [-10, -5, 1, 5, 10]; 
frequencies = [16.77, 8.25, 2.832, 8.12, 16.13]; 

figure
plot(U_out, frequencies); 

% Ukol b
omega_2 = (2.832 * 1e3 / 2500) * 2; 

% Ukol c




