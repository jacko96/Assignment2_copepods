%% Clear
clear all;clc;
close all

%% Parameters - somewhat based on parameter sensitivity analysis
%System parameters
depth=100; %depth of water column [m]
param.dz=1/4; %length or depth of each cell
param.n=depth/param.dz; %number of gridcells
param.z=0.5*param.dz:param.dz:depth-0.5*param.dz; %cell centers (FVM)
time_duration=24*365*1; %[hours] 

%Light parameters
param.I0=200; %[W m^-2] 1 W m^-2 approx. 1 μmole. m2/s ≈ 0.21739130434 W/m2.

%NPD parameters
param.D=(5*(60*60))/10000; %diffusivity [m^2 h^-1]
param.k_w=0.0375; %light absorbtion done by phytoplankton [m^-1]
param.k_c=0.05/1000*14; %selfshading done by phytoplankton and detritus [m^2 (mg N)^-1] [m^2 (mg N)^-1]
param.gmax=0.5/24; % maximum growth rate[h^-1]
param.alpha=0.1/24; % Light sensitivity[m^2 W^-1 h^-1]
param.k_N=0.3*14; % nutrient half-saturation [mg N m^-3]
param.eps=0.03/24; % natural mortality[h^-1]
%param.m=1.5/24*14; % grazing [m^3 (mg N)^-1 h^-1]
param.thao=0.1/24;%0.01; % remineralization[h^-1]
%param.Av=5*10^(-5)*(60*60); % %Vertical diffusion [m^2 h^-1]
param.w=15/24; % sinking of detritus [m h^-1]
param.u=0.001/24; %vertical velocity of phytoplankton [m h^-1]
param.Nb=30; % 30*14 -->bottom concentration of nutrients [mg N m^-3]


%Copepods Parameters
param.S=5; %Number of lifestages for Copepods

param.CtoN=1/5.6;
param.v=0.0052/(24*1000)*param.CtoN; % 0.0052/(24*1000)*param.CtoN-->passive clearance rate coefficient [L mg N^{-3/4} h^{-1}]
param.q=-1/4; %-1/4 -->passive clearance rate exponent (no exponent)
param.h=0.4/(24*1000)*param.CtoN; % maximum ingestion rate coefficent [mg N^{1/4} h^{-1}]
param.h_n=-1/4; %maximum ingestion rate exponent
param.kappa=0.048/(24*1000)*param.CtoN; % respiration rate coefficent [mg N^{1/4} h^{-1}]
param.p=-1/4; %respiration rate exponent
param.assi=0.67; %0.67 -> assimilation rate
param.reproc=((47+0.5)/2)/24; %0.25 -> reproduction efficiency rate
%param.mu_ht=0.003*param.CtoN*0.001; %Mortality by higher trophic levels [mg N^{1/4} mgN^{-2} L^{-2} h^{-1}]
param.Cs=36; %0.001 cm -->Copepod speed velocity [m/h] 
%param.b=((0.5+47)/2)/24/2; %0.37/24 -> reproduction rate [h^{-1}]
%param.mu=1/(((14+365)/2)/24); %1/(((14+365)/2)/24) -> death or mortality rate of copepods [h^-1] --> based on averages found online
param.mu=0.003*param.CtoN/24;

%Copepods grid space
param.m0=0.01; % mgN pr individual
param.ma=1; % mgN pr individual
param.m_bound=exp(linspace(log(param.m0),log(param.ma),param.S+1)); % mgN pr individual
param.m_center=exp(((log(param.m_bound(2:param.S+1))+log(param.m_bound(1:param.S)))/2)); % mgN pr individual
%param.m_ratio=param.m_center(1)./param.m_center(2);
param.m_ratio=param.m_bound(1)./param.m_bound(2);
param.ds=log(param.m_center(2:5))-log(param.m_center(1:4));
param.half_ds=log(param.m_center(1))-log(param.m_bound(1));


%%  ------------------------------ Initialization and ODE-solver  ------------------------------ 

% --------------------------------Initial conditions------------------------------------
C1=0.*ones(1,param.n);
C2=0.*ones(1,param.n);
C3=0.*ones(1,param.n);
C4=0.*ones(1,param.n);
C5=0.*ones(1,param.n);

C1(1:20)=5;
C2(1:20)=5;
C3(1:20)=5;
C4(1:20)=5;
C5(1:20)=5;

P =0*ones(1,param.n);
P(1:floor(param.n/4))=10;
N = 0*ones(1,param.n);
N(1:5)=0.01;
N(end-5:end)=0.01;
D = 5.*ones(1,param.n); %To equal the abundance of copepods pr group

%N(1:end)=0+param.Nb/length(param.z)*param.z;

Initial=[P N D C1 C2 C3 C4 C5];

%  ------------------------------------ ODE ------------------------------------
[t,y]=ode45(@deriv,0:time_duration,Initial,[],param);

% ------------------------------Assign state variables------------------
P = y(:,1:param.n);
N = y(:,param.n+1:2*param.n);
D = y(:,2*param.n+1:3*param.n);

C1 =y(:,3*param.n+1:4*param.n);
C2 =y(:,4*param.n+1:5*param.n);
C3 =y(:,5*param.n+1:6*param.n);
C4 =y(:,6*param.n+1:7*param.n);
C5 =y(:,7*param.n+1:8*param.n);



%% Light and limiting conditions
I=calclight(t,P(end,:),D(end,:),param);

sigma_N=(N(end,:)./(param.k_N+N(end,:)));
sigma_L=(param.alpha.*I(:,end))./(sqrt(param.gmax^2+param.alpha^2.*I(:,end).^2));
%sigma_dot=(sigma_L'.*sigma_N);
sigma=min(sigma_L',sigma_N);

figure
plot(sigma_N,-param.z','Color','k',LineWidth=2)
hold on
plot(sigma_L,-param.z','Color','b',LineWidth=2)
hold on
plot(sigma,-param.z','r*')
legend('Nutrient limiting','Light limiting','Limiting factor','Location','best')
ylabel('Depth [m]')
xlabel('Limiting factor')
grid on; grid minor
box on


%%  ------------------------------ NPD Dynamics  ------------------------------ 

figure('Name', "NPZD Dynamics");

% Subplot for Phytoplankton (P)
subplot(3, 3, 1);
surface(t/24/365, -param.z', P');
shading interp;
title('Phytoplankton');
colorbar;
grid on; grid minor;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);
%clim([0 2])

% Subplot for Nutrients (N)
subplot(3, 3, 2);
surface(t/24/365, -param.z', N');
shading interp;
title('Nutrients');
colorbar;
grid on; grid minor;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);

% Subplot for Detritus (D)
subplot(3, 3, 3);
surface(t/24/365, -param.z', D');
shading interp;
title('Detritus');
colorbar;
grid on; grid minor;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);
%clim([0 100])

% Subplot for C1
subplot(3, 3, 4);
surface(t/24/365, -param.z', C1');
shading interp;
title('C1');
colorbar;
grid on; grid minor;
ylabel('Depth [m]');
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);
%clim([1e-9 1e-4])
ylim([-6 0])

% Subplot for C2
subplot(3, 3, 5);
surface(t/24/365, -param.z', C2');
shading interp;
title('C2');
colorbar;
grid on; grid minor;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);
%clim([1e-9 1e-4])
ylim([-6 0])

% Subplot for C3
subplot(3, 3, 6);
surface(t/24/365, -param.z', C3');
shading interp;
title('C3');
colorbar;
grid on; grid minor;
xlabel('Time [yr]');
%clim([1e-9 1e-4])
ylim([-6 0])

% Subplot for C4
subplot(3, 3, 7);
surface(t/24/365, -param.z', C4');
shading interp;
title('C4');
colorbar;
grid on; grid minor;
xlabel('Time [yr]');
%clim([1e-9 1e-4])
ylim([-6 0])

% Subplot for C5
subplot(3, 3, 8);
surface(t/24/365, -param.z', C5');
shading interp;
title('C5');
colorbar;
grid on; grid minor;
xlabel('Time [yr]');
%clim([1e-9 1e-4])
ylim([-6 0])

% Adjust the layout
sgtitle('State Variables Dynamics [mgN]');
linkaxes(findall(gcf, 'type', 'axes'), 'x');


%%  ------------------------------ Last time step  ------------------------------ 
figure('Name', "NPZD Dynamics");

% Subplot for Phytoplankton (P)
subplot(3, 3, 1);
plot(P(end,:), -param.z','Color','b',LineWidth=2);
title('Phytoplankton');
grid on; grid minor;
%xlabel('Abundance [#]');
%ylim([-10 0])

% Subplot for Nutrients (N)
subplot(3, 3, 2);
plot(N(end,:), -param.z','Color','b',LineWidth=2);
title('Nutrients');
grid on; grid minor;
%xlabel('Abundance [#]');

% Subplot for Detritus (D)
subplot(3, 3, 3);
plot(D(end,:), -param.z','Color','b',LineWidth=2);
title('Detritus');
grid on; grid minor;
%xlabel('Abundance [#]');

% Subplot for C1
subplot(3, 3, 4);
plot(C1(end,:), -param.z','Color','b',LineWidth=2);
title('C1');
grid on; grid minor;
ylabel('Depth [m]');
%xlabel('Abundance [#]');
ylim([-10 0])

% Subplot for C2
subplot(3, 3, 5);
plot(C2(end,:), -param.z','Color','b',LineWidth=2);
title('C2');
grid on; grid minor;
%xlabel('Abundance [#]');
ylim([-10 0])

% Subplot for C3
subplot(3, 3, 6);
plot(C3(end,:), -param.z','Color','b',LineWidth=2);
title('C3');
grid on; grid minor;
xlabel('Abundance [mg N]');
ylim([-10 0])


% Subplot for C4
subplot(3, 3, 7);
plot(C4(end,:), -param.z','Color','b',LineWidth=2);
title('C4');
grid on; grid minor;
xlabel('Abundance [mg N]');
ylim([-10 0])

% Subplot for C5
subplot(3, 3, 8);
plot(C5(end,:), -param.z','Color','b',LineWidth=2);
title('C5');
grid on; grid minor;
xlabel('Abundance [mg N]');
ylim([-10 0])


% Adjust the layout
% sgtitle('State Variables Dynamics [mgN]');
% linkaxes(findall(gcf, 'type', 'axes'), 'x');

%% Somatic growth
Set_depth = round(20); %Sets the depth at which you want to investigate somatic growth
Cop_space = [C1(:,Set_depth) C2(:,Set_depth) C3(:,Set_depth) C4(:,Set_depth) C5(:,Set_depth)]; %bottom of water column
figure()
h = surface(t/24/365,log10(param.m_center),Cop_space');
set(gca, 'ZScale', 'log'); % Set z-axis to logarithmic scale
shading interp
colorbar
title(['Somatic growth at ' num2str(Set_depth) 'm depth']);
xlabel('Time [yr]')
ylabel('Body mass [mgN/#]')

% Plot ylines and add red text "group x" beside each yline
yline(log10(param.m_center(1)),'-.k','LineWidth',2)
text(0.5, param.m_center(1), 'Group 1', 'HorizontalAlignment', 'right', 'Color', 'k');

yline(log10(param.m_center(2)),'-.k','LineWidth',2)
text(0.5, param.m_center(2), 'Group 2', 'HorizontalAlignment', 'left', 'Color', 'k');

yline(log10(param.m_center(3)),'-.k','LineWidth',2)
text(0.5, param.m_center(3), 'Group 3', 'HorizontalAlignment', 'left', 'Color', 'k');

yline(log10(param.m_center(4)),'-.k','LineWidth',2)
text(0.5, param.m_center(4), 'Group 4', 'HorizontalAlignment', 'left', 'Color', 'k');

yline(log10(param.m_center(5)),'-.k','LineWidth',2)
text(0.5, param.m_center(5), 'Group 5', 'HorizontalAlignment', 'left', 'Color', 'k');


%%  ----------------------------- Grid sensitivity analysis 
figure('Name',"Grid Sensitivity Analysis");
param.dz=1;

% Define different values of param.dz
deltaz_values = [param.dz/4,param.dz/2,param.dz,param.dz*2,param.dz*4];

for j = 1:length(deltaz_values)
    % Update param.dz
    param.dz = deltaz_values(j);

    % Recalculate the number of grid cells and cell interval
    param.n = floor(depth / param.dz);
    param.z = 0.5 * param.dz : param.dz : depth - 0.5 * param.dz;

    % Initialize arrays to store results
    C1=0.*ones(1,param.n);
    C2=0.*ones(1,param.n);
    C3=0.*ones(1,param.n);
    C4=0.*ones(1,param.n);
    C5=0.*ones(1,param.n);

    C1(1:20)=5;
    C2(1:20)=5;
    C3(1:20)=5;
    C4(1:20)=5;
    C5(1:20)=5;

    P =0*ones(1,param.n);
    P(1:floor(param.n/4))=10;
    N = 0*ones(1,param.n);
    N(1:5)=0.01;
    N(end-5:end)=0.01;
    D = 5.*ones(1,param.n); %To equal the abundance of copepods pr group

    Initial2=[P(1,:) N(1,:) D(1,:) C1(1,:) C2(1,:) C3(1,:) C4(1,:) C5(1,:)];

    % Solve ODEs
    [t2, y2] = ode45(@deriv, 0:time_duration, Initial2,[], param);

    % Extract P, N, D for the current grid spacing
    P2 = y2(:,1:param.n);
    N2 = y2(:,param.n+1:2*param.n);
    D2 = y2(:,2*param.n+1:3*param.n);

    C1 =y2(:,3*param.n+1:4*param.n);
    C2 =y2(:,4*param.n+1:5*param.n);
    C3 =y2(:,5*param.n+1:6*param.n);
    C4 =y2(:,6*param.n+1:7*param.n);
    C5 =y2(:,7*param.n+1:8*param.n);

    % Plot the results with varying color intensity
    color_index = (j - 1) / (length(deltaz_values) - 1); % Gradient index
    color = [1-color_index, color_index, 0]; % RGB color based on the index
    plot(C5(end, :), -param.z, 'LineWidth', 2, 'Color', color, 'DisplayName', ['\Delta z = ', num2str(param.dz)]);
    drawnow
    hold on;
end

title('Grid Sensitivity Analysis');
xlabel('Copepods group 5 abundance [mg N/m^3]');
ylabel('Depth [m]');
legend('Location', 'best');
grid on; grid minor;




%% ------------------------ Parameter sensitivity analysis (and somatic growths)

depth = 100; % Depth of water column [m]
param.dz = 1/4; % Length or depth of each cell
param.n = depth / param.dz; % Number of grid cells
param.z = 0.5 * param.dz:param.dz:depth - 0.5 * param.dz; % Cell interval


Set_depth = round(1/param.dz); %Sets the depth at which you want to investigate somatic growth


% Define different values of specific parameter
val= [param.k_N*10 param.k_N*2, param.k_N, param.k_N/2, param.k_N/10]; %Change/update the parameter inside this vector
% use either param.k_N, param.v, param.Cs, param.h, or param.eps, 
% Define colormap
colormap('default')

%figure('Name', 'Sensitivity analysis')

% Iterate over each value of specific parameter
for i = 1:length(val)
    % Update parameter value in param struct
    param.k_N = val(i); %Change/update this parameter

    % Initialize arrays to store results
    C1=0.*ones(1,param.n);
    C2=0.*ones(1,param.n);
    C3=0.*ones(1,param.n);
    C4=0.*ones(1,param.n);
    C5=0.*ones(1,param.n);

    C1(1:20)=5;
    C2(1:20)=5;
    C3(1:20)=5;
    C4(1:20)=5;
    C5(1:20)=5;

    P =0*ones(1,param.n);
    P(1:floor(param.n/4))=10;
    N = 5*ones(1,param.n);
    N(1:5)=param.Nb;
    N(end-5:end)=param.Nb;
    D = 5.*ones(1,param.n); %To equal the abundance of copepods pr group

    Initial3=[P(1,:) N(1,:) D(1,:) C1(1,:) C2(1,:) C3(1,:) C4(1,:) C5(1,:)];

    % Call ode23 solver with the correct initial conditions
    [t, y1] = ode45(@deriv, 0:time_duration, Initial3, [], param);

    % Extract P, N, D for the current value of gmax
    P1 = y1(:, 1:param.n);
    N1 = y1(:, param.n + 1:2 * param.n);
    D1 = y1(:, 2 * param.n + 1:3 * param.n);
    C1 = y1(:, 3*param.n+1:4*param.n);
    C2 = y1(:, 4*param.n+1:5*param.n);
    C3 = y1(:, 5*param.n+1:6*param.n);
    C4 = y1(:, 6*param.n+1:7*param.n);
    C5 = y1(:, 7*param.n+1:8*param.n);

    % Define color based on gradient index
    color_index = (i - 1) / (length(val) - 1); % Gradient index
    color = [1-color_index, color_index, 0]; % RGB color based on the index

    % Plot P at the end of the simulation with gradient color
    figure(1)
    plot(C5(end, :), -param.z, 'MarkerSize', 15, 'LineWidth', 2.5, 'Color', color, 'LineStyle', '-', 'DisplayName', "k_{N}"+" = " + num2str(val(i)) + " mg N m^{-3} ");
    %hold on
    %plot(P(end, :), -param.z, 'MarkerSize', 15, 'LineWidth', 2.5, 'Color', color, 'LineStyle', '--', 'DisplayName', "C1" + " = " + num2str(val(i)) + " ");
    drawnow
    hold on;
    legend('Location', 'southeast')
    grid on; grid minor
    xlabel('# of copepods')
    ylabel('Depth [m]')
    
    figure(2)
    subplot(3,2,i)
        switch i
            case 1
                Cop_space = [C1(:,Set_depth) C2(:,Set_depth) C3(:,Set_depth) C4(:,Set_depth) C5(:,Set_depth)]; %bottom of water column
                h = surface(t/24/365,log10(param.m_center),Cop_space');
                set(gca, 'ZScale', 'log'); % Set z-axis to logarithmic scale
                shading interp
                colorbar
                title(['Somatic growth at ' num2str(Set_depth) 'm depth for k_{N} =' num2str(val(i)) ' mg N m^{-3}']);
                xlabel('Time [yr]')
                ylabel('Body mass [mgN/#]')
                % Plot ylines and add red text "group x" beside each yline
                yline(log10(param.m_center(1)),'-.k','LineWidth',1)
                text(0.5, param.m_center(1), 'Group 1', 'HorizontalAlignment', 'right', 'Color', 'k');
                yline(log10(param.m_center(2)),'-.k','LineWidth',1)
                text(0.5, param.m_center(2), 'Group 2', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(3)),'-.k','LineWidth',1)
                text(0.5, param.m_center(3), 'Group 3', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(4)),'-.k','LineWidth',1)
                text(0.5, param.m_center(4), 'Group 4', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(5)),'-.k','LineWidth',1)
                text(0.5, param.m_center(5), 'Group 5', 'HorizontalAlignment', 'left', 'Color', 'k');
            case 2
                Cop_space = [C1(:,Set_depth) C2(:,Set_depth) C3(:,Set_depth) C4(:,Set_depth) C5(:,Set_depth)]; %bottom of water column
                h = surface(t/24/365,log10(param.m_center),Cop_space');
                set(gca, 'ZScale', 'log'); % Set z-axis to logarithmic scale
                shading interp
                colorbar
                title(['Somatic growth at ' num2str(Set_depth) 'm depth for k_{N} =' num2str(val(i)) ' mg N m^{-3}']);
                xlabel('Time [yr]')
                ylabel('Body mass [mgN/#]')
                % Plot ylines and add red text "group x" beside each yline
                yline(log10(param.m_center(1)),'-.k','LineWidth',1)
                text(0.5, param.m_center(1), 'Group 1', 'HorizontalAlignment', 'right', 'Color', 'k');
                yline(log10(param.m_center(2)),'-.k','LineWidth',1)
                text(0.5, param.m_center(2), 'Group 2', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(3)),'-.k','LineWidth',1)
                text(0.5, param.m_center(3), 'Group 3', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(4)),'-.k','LineWidth',1)
                text(0.5, param.m_center(4), 'Group 4', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(5)),'-.k','LineWidth',1)
                text(0.5, param.m_center(5), 'Group 5', 'HorizontalAlignment', 'left', 'Color', 'k');
            case 3
                Cop_space = [C1(:,Set_depth) C2(:,Set_depth) C3(:,Set_depth) C4(:,Set_depth) C5(:,Set_depth)]; %bottom of water column
                h = surface(t/24/365,log10(param.m_center),Cop_space');
                set(gca, 'ZScale', 'log'); % Set z-axis to logarithmic scale
                shading interp
                colorbar
                title(['Somatic growth at ' num2str(Set_depth) 'm depth for k_{N} =' num2str(val(i)) ' mg N m^{-3}']);
                xlabel('Time [yr]')
                ylabel('Body mass [mgN/#]')
                % Plot ylines and add red text "group x" beside each yline
                yline(log10(param.m_center(1)),'-.k','LineWidth',1)
                text(0.5, param.m_center(1), 'Group 1', 'HorizontalAlignment', 'right', 'Color', 'k');
                yline(log10(param.m_center(2)),'-.k','LineWidth',1)
                text(0.5, param.m_center(2), 'Group 2', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(3)),'-.k','LineWidth',1)
                text(0.5, param.m_center(3), 'Group 3', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(4)),'-.k','LineWidth',1)
                text(0.5, param.m_center(4), 'Group 4', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(5)),'-.k','LineWidth',1)
                text(0.5, param.m_center(5), 'Group 5', 'HorizontalAlignment', 'left', 'Color', 'k');
            case 4
                Cop_space = [C1(:,Set_depth) C2(:,Set_depth) C3(:,Set_depth) C4(:,Set_depth) C5(:,Set_depth)]; %bottom of water column
                h = surface(t/24/365,log10(param.m_center),Cop_space');
                set(gca, 'ZScale', 'log'); % Set z-axis to logarithmic scale
                shading interp
                colorbar
                title(['Somatic growth at ' num2str(Set_depth) 'm depth for k_{N} =' num2str(val(i)) ' mg N m^{-3}']);
                xlabel('Time [yr]')
                ylabel('Body mass [mgN/#]')
                % Plot ylines and add red text "group x" beside each yline
                yline(log10(param.m_center(1)),'-.k','LineWidth',1)
                text(0.5, param.m_center(1), 'Group 1', 'HorizontalAlignment', 'right', 'Color', 'k');
                yline(log10(param.m_center(2)),'-.k','LineWidth',1)
                text(0.5, param.m_center(2), 'Group 2', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(3)),'-.k','LineWidth',1)
                text(0.5, param.m_center(3), 'Group 3', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(4)),'-.k','LineWidth',1)
                text(0.5, param.m_center(4), 'Group 4', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(5)),'-.k','LineWidth',1)
                text(0.5, param.m_center(5), 'Group 5', 'HorizontalAlignment', 'left', 'Color', 'k');
            case 5
                Cop_space = [C1(:,Set_depth) C2(:,Set_depth) C3(:,Set_depth) C4(:,Set_depth) C5(:,Set_depth)]; %bottom of water column
                h = surface(t/24/365,log10(param.m_center),Cop_space');
                set(gca, 'ZScale', 'log'); % Set z-axis to logarithmic scale
                shading interp
                colorbar
                title(['Somatic growth at ' num2str(Set_depth) 'm depth for k_{N} =' num2str(val(i)) ' mg N m^{-3}']);
                xlabel('Time [yr]')
                ylabel('Body mass [mgN/#]')
                % Plot ylines and add red text "group x" beside each yline
                yline(log10(param.m_center(1)),'-.k','LineWidth',1)
                text(0.5, param.m_center(1), 'Group 1', 'HorizontalAlignment', 'right', 'Color', 'k');
                yline(log10(param.m_center(2)),'-.k','LineWidth',1)
                text(0.5, param.m_center(2), 'Group 2', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(3)),'-.k','LineWidth',1)
                text(0.5, param.m_center(3), 'Group 3', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(4)),'-.k','LineWidth',1)
                text(0.5, param.m_center(4), 'Group 4', 'HorizontalAlignment', 'left', 'Color', 'k');
                yline(log10(param.m_center(5)),'-.k','LineWidth',1)
                text(0.5, param.m_center(5), 'Group 5', 'HorizontalAlignment', 'left', 'Color', 'k');
        end


end



function I=calclight(t,P,D,param)
light=param.I0-param.I0.*sin(2*pi/(24*365).*t);
integral=cumsum((param.k_w+param.k_c.*(P+D)).*param.dz+1/2*param.dz.*(param.k_w+param.k_c.*(P+D)));
%I=light'.*exp(-param.k_w.*param.z-param.k_c.*integral'); %I=param.I0*exp(-integral');
I=light'.*exp(-integral');
end

function dydt = deriv(t,y,param) %calculates the dynamics of the inputs throughout the water column

P = y(1:param.n);
N = y(param.n+1:2*param.n);
D = y(2*param.n+1:3*param.n);

C1 = y(3*param.n+1:4*param.n);
C2 = y(4*param.n+1:5*param.n);
C3 = y(5*param.n+1:6*param.n);
C4 = y(6*param.n+1:7*param.n);
C5 = y(7*param.n+1:8*param.n);

%-------------------- Setting negative values to zero
% P(P<0)=0; N(N<0)=0; D(D<0)=0; C1(C1<0)=0; C2(C2<0)=0; C3(C3<0)=0;
% C4(C4<0)=0; C5(C5<0)=0;

% ---------------------------- Fluxes

ix=2:param.n; %Vector notation for UDS
ixx=2:param.n-1; %Vector notation for CDS

%Phytoplankton
JaP(ix)=0;
JdP(ix)=-param.D*(P(ix)-P(ix-1))/param.dz;

%Nutrients
JaN(ix)=0;
JdN(ix)=-param.D*(N(ix)-N(ix-1))/param.dz;

%Detritus
JaD(ix)=-param.w*(D(ix-1)); 
JdD(ix)=-param.D*(D(ix)-D(ix-1))/param.dz;

%Copepods
JaC(ix)=0;
%JdC(ix)=param.Cs*(P(ix)-P(ix-1))/param.dz;
JdC(ixx) = param.Cs * (P(ixx+1) - P(ixx-1)) / (2 * param.dz); % Using CDS to get the Copepods to move everywhere in the water column

% --------------------------------- Boundary conditions and alterations
%Phytoplankton
JaP(1)=0;
JaP(param.n+1)=0;
JdP(1)=0; % or -param.D*(P(1)-0)/param.dz;
JdP(param.n+1)=0;

%Nutrients
JaN(1)=0;
JaN(param.n+1)=0; 
JdN(1)=0;
JdN(param.n+1)=-param.D*(param.Nb-N(end))/param.dz; % eller -param.D*(param.Nb-N(param.n))/param.dz; 

%Detritus
JaD(1)=0;
JaD(param.n+1)=-param.w*D(end); 
JdD(1)=0;
JdD(param.n+1)=0;

%Copepods
JaC(1)=0;
JaC(param.n+1)=0;
JdC(1)=0;%param.Cs*(P(2)-P(1))/(param.dz);% 0 ikke sikker på, om det er helt rigtigt
JdC(param.n+1)=0;


%JP=JaP+JdP;
%JN=JaN+JdN;
%JD=JaD+JdD;
%JC=JaC+JdC;

%Assembling the equations
I = calclight(t,P,D,param);

sigma_N=(N./(param.k_N+N));
sigma_L=(param.alpha.*I)./(sqrt(param.gmax^2+param.alpha^2.*I.^2));
%sigma=(sigma_L'.*sigma_N);
sigma=min(sigma_L',sigma_N);

% ---------------------------- Advection-Diffusion
gradientP=(-(JaP(2:param.n+1)+JdP(2:param.n+1)-JaP(1:param.n)-JdP(1:param.n))/param.dz)';
gradientN=(-(JaN(2:param.n+1)+JdN(2:param.n+1)-JaN(1:param.n)-JdN(1:param.n))/param.dz)';
gradientD=(-(JaD(2:param.n+1)+JdD(2:param.n+1)-JaD(1:param.n)-JdD(1:param.n))/param.dz)';
gradientC=(-(JaC(2:param.n+1)+JdC(2:param.n+1)-JaC(1:param.n)-JdC(1:param.n))/param.dz)';

% ----------------------------- Copepods

% Feeding level (0-1):
f = (param.v.*(param.m_center.^param.q).*P)./(param.v.*(param.m_center.^param.q).*P+param.h.*param.m_center.^param.h_n);
%f =(param.v.*(param.m_center.^param.q).*P)./(param.v.*(param.m_center.^param.q).*P + param.h.*param.m_center.^param.h_n); % P/1000 fordi der er 1000 L på en m3.

fc = param.kappa./(param.assi*param.h).*param.m_center.^(param.p - param.h_n);
% biomass accumulation/production rate:
v = param.assi*param.h.*param.m_center.^param.h_n.*(f - fc);
g = max(0,v); % net energy gain
mu_st = min(0,v); %mu constant?
mu = mu_st-param.mu; % copepod mortality
b = param.reproc.*g(:,end); % copepod birth rate
gamma = (g - mu)./(1 - (param.m_ratio).^(1 - mu./g)); % Overgang fra lille til stor copepod
gamma(gamma<0)=0;
%gamma = max(zeros(size(gamma)),gamma); % hvis gamma bliver negativ

% dC1dt = b.*C5 + g(:,1).*C1 - gamma(:,1).*C1 - mu(:,1).*C1 + gradientC;
% dC2dt = gamma(:,1).*C1 + g(:,2).*C2 - gamma(:,2).*C2 - mu(:,2).*C2 + gradientC;
% dC3dt = gamma(:,2).*C2 + g(:,3).*C3 - gamma(:,3).*C3 - mu(:,3).*C3 + gradientC;
% dC4dt = gamma(:,3).*C3 + g(:,4).*C4 - gamma(:,4).*C4 - mu(:,4).*C4 + gradientC;
% dC5dt = gamma(:,4).*C4 - mu(:,5).*C5 + gradientC;

% Compute the rates of change dC1dt to dC5dt
dC1dt = b.*C5 + g(:,1).*C1 - gamma(:,1).*C1 - mu(:,1).*C1 + gradientC.*(C1 > 0);
dC2dt = gamma(:,1).*C1 + g(:,2).*C2 - gamma(:,2).*C2 - mu(:,2).*C2 + gradientC.*(C2 > 0);
dC3dt = gamma(:,2).*C2 + g(:,3).*C3 - gamma(:,3).*C3 - mu(:,3).*C3 + gradientC.*(C3 > 0);
dC4dt = gamma(:,3).*C3 + g(:,4).*C4 - gamma(:,4).*C4 - mu(:,4).*C4 + gradientC.*(C4 > 0);
dC5dt = gamma(:,4).*C4 - mu(:,5).*C5 + gradientC.*(C5 > 0);

Copepod_loss=mu(:,1).*C1+mu(:,2).*C2+mu(:,3).*C3+mu(:,4).*C4+mu(:,5).*C5;


% ----------------------------- NPD equations
Phytoplankton_loss=param.h.*param.m_center.^param.h_n.*(f);
Phytoplankton_loss2=Phytoplankton_loss(:,1).*C1./param.m_center(1)+Phytoplankton_loss(:,2).*C2./param.m_center(2)+Phytoplankton_loss(:,3).*C3./param.m_center(3)+Phytoplankton_loss(:,4).*C4./param.m_center(4)+Phytoplankton_loss(:,5).*C5./param.m_center(5);

dPdt = param.gmax .* sigma .* P - param.eps .* P -Phytoplankton_loss2+ gradientP; %new equation utilizing equations from the paper
%dPdt = param.gmax .* sigma .* P - param.eps .* P - sum(param.v.*(param.m_center.^param.q).*P/1000,2) + gradientP; % Har erstattet grazing rate med encounter rate fra ligning (1)
dNdt = -param.gmax .* sigma .* P + param.thao .* D + gradientN; %new equation utilizing equations from the paper
dDdt = param.eps .* P  - param.thao .* D+Copepod_loss+gradientD; %new equation utilizing equations from the paper
%dDdt = param.eps .* P  - param.thao.*D + param.m.*(P.*P) +param.thao*Phytoplankton_loss2+gradientD; %new equation utilizing equations from the paper

%adding the loss of copepods to detritus equation



dydt=[dPdt;dNdt;dDdt;dC1dt;dC2dt;dC3dt;dC4dt;dC5dt];

end