%% Parameters for arap baseline.
clc; clear all; close all;
model_name = 'square_21_spikes';
algorithm = 'arap';
iter = 100;
% Run terminal command.
command = ['./build/demo_bin ', './model/', model_name, ' ', algorithm, ' ', num2str(iter)];
system(command);

% Create data file names.
clear arap_*;
file_name = [algorithm, '-', num2str(iter), '.txt'];
arap_data = readtable(['./data/', model_name, '/', file_name], 'Delimiter', '\t');
arap_headers = arap_data.Properties.VariableNames;
for i = 1 : length(arap_headers)
    eval(['arap_', lower(arap_headers{i}), '=arap_data{:, ', num2str(i), '};']);
end

%% Test different algorithms.
algorithm = 'admm-fixed';
rho = 0.00087;
command = ['./build/demo_bin ', './model/', model_name, ' ', algorithm, ' ', num2str(iter), ' ', num2str(rho)];
system(command);

% Get data files.
clear admm_*;
file_name = [algorithm, '-', num2str(iter), '-', num2str(rho), '.txt'];
admm_data = readtable(['./data/', model_name, '/', file_name], 'Delimiter', '\t');
admm_headers = admm_data.Properties.VariableNames;
for i = 1 : length(admm_headers)
    eval(['admm_', lower(admm_headers{i}), '=admm_data{:, ', num2str(i), '};']);
end

%% Plot.
subplot(3, 1, 1);
%conv_arap = (arap_total(2:end)./arap_total(1:end-1));
%conv_admm = (admm_arap(2:end)./admm_arap(1:end-1));
%plot(arap_iteration(1:end-1), log(conv_arap), 'r', admm_iteration(1:end-1), log(conv_admm), 'b');
plot(arap_iteration, arap_total, 'r', admm_iteration, admm_arap, 'b');
legend('arap energy', 'admm energy');
xlabel('iteration');
ylabel('arap/admm energy');

subplot(3, 1, 2);
plot(arap_iteration, arap_xnorm, 'r', admm_iteration, admm_xnorm, 'b');
legend('arap norm', 'admm norm');
xlabel('iteration');
ylabel('solution norm');

subplot(3, 1, 3);
plot(arap_iteration, arap_xdiffnorm, 'r', admm_iteration, admm_xdiffnorm, 'b');
legend('arap xdiff norm', 'admm xdiff norm');
xlabel('iteration');
ylabel('diff solution norm');

% Analyse convergence rates.
figure;
subplot(3, 1, 1);
plot(admm_iteration, admm_rho, 'b');
xlabel('iteraiton');
ylabel('rho');

subplot(3, 1, 2);
plot(admm_iteration, admm_rotationavg, 'b');
xlabel('iteration');
ylabel('rotation error');

if exist('admm_vertexavg', 'var')
    subplot(3, 1, 3);
    plot(admm_iteration, admm_vertexavg, 'b');
    xlabel('iteraiton');
    ylabel('vertex error');
    clear admm_vertex;
end