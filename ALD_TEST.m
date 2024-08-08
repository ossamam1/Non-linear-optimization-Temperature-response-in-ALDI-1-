clear *; close all; clc

%% Nonlinear optimization for the reactor furnace set temperature via CasADI


load ALD_vacuum_data.mat
input_data = [data{:,[1,2]} data{:,1}./data{:,2}];
output_data = data{:,end};
sys = getenv('USERNAME');
path = fullfile('\\home.org.aalto.fi\ossamam1\data\Desktop\casADi');
addpath(genpath(path));
theta = size(input_data,2);
import casadi.* 
opti = casadi.Opti();
W1 = opti.variable(theta);
b1 = opti.variable();
y_calc = (input_data*W1+b1); 
obj = norm(y_calc - output_data);
opti.minimize(obj^2) 
opti.set_initial(W1,2*rand(theta,1)-1)
opti.solver('ipopt')
opti.solve()
W1_sol = opti.value(W1); 
b1_sol = opti.value(b1);
y_pred = input_data*W1_sol+b1_sol;
h1 = plot((output_data),(y_pred),'ko');
xlabel('Experimental reactor furnace set point (°C)')
ylabel('Model predicted reactor furnace set point (°C)')
title ('Experimental vs model predicted values')
xlim([0 max(output_data)]) 
refline
figure()
plot(output_data,output_data-y_pred,'ko')
xlabel('Experimental reactor furnace set point (°C)')
ylabel('Experimental - Predicted values')
xlim([0 max(output_data)])
n_train=size(input_data,1);
rmse = sqrt( sum( (output_data - y_pred).^2 ) / n_train ) 
mean_error = mae(output_data,y_pred)
%% Mesh plots
figure()
[X,Y]=meshgrid(input_data(:,1),input_data(:,2));
z = W1_sol(1).*X + W1_sol(2).*Y + W1_sol(3).*X./Y + b1_sol;
surf(X,Y,z)
c = colorbar;
c.Label.String = 'Temperature (°C)';
c.Label.FontSize = 12;
colormap(turbo);
xlabel('Connection hub set point (°C)')
ylabel('Average bed temperature (°C)')
zlabel('Reactor furnace set point (°C)') 
title ('Non linear optimization via CasADI')
xlim([0 200])
ylim([0 800])
zlim([0 800])