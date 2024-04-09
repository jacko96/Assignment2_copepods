% Exercise 4
clear;close all;clc;

% Define parameters
run('Parameters.m')

% Defining grid
param.n = 100; % no. cells
param.dz = param.D/param.n; % grid size
param.z =  param.dz/2:param.dz:param.D; % grid

% 5. solve
P0 = ones(param.n,1)*10; % Initial conditions:
N0 = ones(param.n,1)*5;
D0 = ones(param.n,1);
Y0 = [P0;N0;D0];
tspan = 0:1:2500;
[t,Y] = ode45(@(t,Y)odefun4(t,Y,param),tspan,Y0);
[~, pI,pN] = odefun4(t(end),Y(end,:),param); % Finding the limiting factors (light or minerals)

Y = Y';
P = Y(1:param.n,:);
N = Y(param.n+1:(2*param.n),:);
D = Y((2*param.n+1):end,:);


%% 6. plot
close all
figure()
subplot(3,1,1)
% Plankton
%figure(1)
image(P,'CDataMapping','scaled')
colorbar
xlabel('Time [100 hours]')
ylabel('Depth [meters]')
title('Phytoplankton conc. [cells/m3]')

subplot(3,1,2)
% Nutrients
%figure(2)
image(N,'CDataMapping','scaled')
colorbar
xlabel('Time [hours]')
ylabel('Depth [meters]')
title('Nutrient conc. [mmol nutrient/m3]')

subplot(3,1,3)
% Detritus
%figure(3)
image(D,'CDataMapping','scaled')
colorbar
xlabel('Time [hours]')
ylabel('Depth [meters]')
title('Detritus conc. [mmol detritus/m3]')


%% Limiting factors
figure(4)
plot(pI,param.z)
set(gca, 'YDir','reverse')
hold on
plot(pN,param.z)
set(gca, 'YDir','reverse')
hold off
legend('Light','Minerals')
title('Growth limiting factors')
ylabel('Depth [meters]')
xlabel('Growth fraction')

%% Sensitivity analysis
close all
% Increasing grazing by 10
param.gamma = 10*0.5/24; % Grazing (m3 (mmol N)^(-1) 1/h)
% Increasing diffusivity by 10
%param.d = 10*3*10^(-4)/(1/(60*60)); % Diffusion constant (m2/h)

[t,Y_10gamma] = ode45(@(t,Y)odefun4(t,Y,param),tspan,Y0);
[~, pI_10gamma,pN_10gamma] = odefun4(t(end),Y_10gamma(end,:),param); % Finding the limiting factors (light or minerals)
Y_10gamma = Y_10gamma';
P_10gamma = Y_10gamma(1:param.n,:);
N_10gamma = Y_10gamma(param.n+1:(2*param.n),:);
D_10gamma = Y_10gamma((2*param.n+1):end,:);

figure()
plot(P(:,end),param.z)
hold on
plot(P_10gamma(:,end),param.z)
xlabel('Phytoplankton conc. [cells/m3]')
ylabel('Depth [meters]')
legend('Original','10x Diffusivity')
%%
[pmax,zmax] = max(P(:,end));
sigmamax = std(P(:,end));

[pmax_new,zmax_new] = max(P_10gamma(:,end));
sigmamax_new = std(P_10gamma(:,end));

