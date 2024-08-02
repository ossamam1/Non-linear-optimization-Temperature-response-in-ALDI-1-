clear *; close all; clc

%% Nonlinear optimization for the reactor furnace set temperature via CasADI

% Read data from excel
% data = readtable("ALDI1-0005.xlsx",'Sheet','Sheet1','Range','I1:K54');
% Save data to mat file
% Acknowledgement to Qasim Riaz (RA in Ville group)
% save ALD data
% Load experimental data

load ALD_vacuum_data.mat % ALDI-1 setup data at vacuum conditions
%load ALD_atmospheric_data.mat % ALDI-1 setup data at atmospheric conditions

% Create input output matrices
input_data = [data{:,[1,2]} data{:,1}./data{:,2}]; %Input variables - 1.Connection hub temperature (T1), 2. Average bed temperature (T2), 3. T1/T2
output_data = data{:,end}; % output_data = table2array(data(:,end))

%addpath(genpath('\\home.org.aalto.fi\ossamam1\data\Desktop\casADi')); 
sys = getenv('USERNAME');
path = fullfile('X:\EQUIPMENT_LOGS\ALDI-1\Analyzed data\Matlab code\casADi');
addpath(genpath(path));



theta = size(input_data,2); % Number of descriptors
% casADi implementation
import casadi.* %Adding casadi to MATLAB
opti = casadi.Opti(); % Creates an Opti Stack structure
% Define the model variables which will be regressed to get optimal solution
W1 = opti.variable(theta); % Define weight matrix
b1 = opti.variable(); %Intercept term
% MODEL OUTPUT
y_calc = (input_data*W1+b1); %input_data(:,1)*W1(1) + input_data(:,2)*W1(2) + input_data(:,3)*W1(3) + b1 
% Objective function to be minimized
obj = norm(y_calc - output_data);
opti.minimize(obj^2) 
% Initialization of weight matrix - if bnot specified then initial values taken as 0.Current initialization defines value between -1 and 1
opti.set_initial(W1,2*rand(theta,1)-1) % Sets an initial solution (default: zeros(n,1))
% Define optimization solver - IPOPT/SQPMETHOD
opti.solver('ipopt') % Chooses a solver (IPOPT, qpOASES, ...)
% Solve optimization problem
opti.solve() % Executes the solver
% Convert from casadi matrices to double 
W1_sol = opti.value(W1); % Retrieves the optimal solution of (W1)
b1_sol = opti.value(b1); % Retrieves the optimal solution of (b1)
% Calculate the predicted value based on opti results
y_pred = input_data*W1_sol+b1_sol;
h1 = plot((output_data),(y_pred),'ko');
xlabel('Experimental reactor furnace set point (°C)') % x axis label
ylabel('Model predicted reactor furnace set point (°C)') % y axis label
title ('Experimental vs model predicted values')
xlim([0 max(output_data)]) % x axis limits from zero to maximum point
% Add a reference line to the plot
refline
% Error plot 
figure()
plot(output_data,output_data-y_pred,'ko')
xlabel('Experimental reactor furnace set point (°C)') % x axis label
ylabel('Experimental - Predicted values') % y axis label
xlim([0 max(output_data)]) % x axis limits from zero to maximum point
n_train=size(input_data,1);
% RMSE 
rmse = sqrt( sum( (output_data - y_pred).^2 ) / n_train )
% Mean absolute error 
mean_error = mae(output_data,y_pred)

%% Mesh plots
figure()
[X,Y]=meshgrid(input_data(:,1),input_data(:,2));
z = W1_sol(1).*X + W1_sol(2).*Y + W1_sol(3).*X./Y + b1_sol;
surf(X,Y,z)
c = colorbar;
c.Label.String = 'Temperature (°C)'; % colorbar axis label
c.Label.FontSize = 12; % colorbar axis label font size
colormap(turbo);
%colormap(slanCM('virdis'))

xlabel('Connection hub set point (°C)') % x axis label
ylabel('Average bed temperature (°C)') % y axis label
zlabel('Reactor furnace set point (°C)') % z axis label
title ('Non linear optimization via CasADI')
xlim([0 200]) % x axis limits
ylim([0 800]) % y axis limits
zlim([0 800]) % z axis limits