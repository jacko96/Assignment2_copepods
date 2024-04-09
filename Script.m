clear all
close all

depth=300; %depth of water column [m]
param.dz=1; %length or depth of each cell
param.n=depth/param.dz; %number of gridcells
%param.z=0.5*param.dz:param.dz:param.n-0.5*param.dz; %cell interval
param.z=0.5*param.dz:param.dz:depth-0.5*param.dz; %cell interval


time_duration=45000; %[hours] simulating now 5,13 years
param.I0=200; %[W m^-2] 1 W m^-2 approx. 1 μmole. m2/s ≈ 0.21739130434 W/m2.


param.D=(5*(60*60))/10000; %diffusivity [m^2 h^-1] 
param.k_w=0.0375; %light absorbtion done by phytoplankton [m^-1]
param.k_c=0.05/1000; %selfshading done by phytoplankton and detritus [m^2 (mmol N)^-1] [m^2 (mmol N)^-1]
param.gmax=0.5/24; % maximum growth rate[h^-1]
param.alpha=0.1/24; % Light sensitivity[m^2 W^-1 h^-1]
param.k_N=0.3; % nutrient half-saturation [mmol N m^-3]
param.eps=0.03/24; % natural mortality[h^-1]
param.m=1.5/24; % grazing [m^3 (mmol N)^-1 h^-1]
param.thao=0.1/24;%0.01; % remineralization[h^-1]
%param.Av=5*10^(-5)*(60*60); % %Vertical diffusion [m^2 h^-1]
param.w=15/24; % sinking of detritus [m h^-1]
param.u=0.001/24; %vertical velocity of phytoplankton [m h^-1]
param.Nb=30; %bottom concentration of nutrients [mmol N m^-3]

%param.HI=50*(60*60*24)/1000; %[mmol photons m^-2 day^-1]
%param.HN=0.02; %[mmol nutrient m^-3]

P = 0 * param.z; %
N = 0 * param.z; %
D = 0 * param.z; %
%N(1,end)=30;
N(1:end)=0+param.Nb/length(param.z)*param.z;
P(1,1)=0.1;
%D(end,1)=0;




[t,y]=ode45(@deriv,[0:time_duration],[P N D],[],param);


P=y(:,1:param.n);
N=y(:,param.n+1:2*param.n);
D=y(:,2*param.n+1:3*param.n);



%% ------------------------------- Limiting factors
I=calclight(t(end),P(end,:),D(end,:),param);

% I_min=I(end,:)./max(I(end,:));
% N_min=N(end,:)./max(N(end,:));
% Mini=min(I_min,N_min);
% sigma=I_min.*N_min;
sigma_N_end=(N(end,:)./(param.k_N+N(end,:)));
sigma_L_end=(param.alpha.*I)./(sqrt(param.gmax^2+param.alpha^2.*I.^2));
sigma_end=(sigma_L_end.*sigma_N_end);
sigma_min=min(sigma_L_end,sigma_N_end);

figure('Name',"Limiting factor")
plot(sigma_L_end, -param.z, '-b', 'LineWidth', 3, 'DisplayName', "Light limitation")
hold on
plot(sigma_N_end, -param.z, '-g', 'LineWidth', 3, 'DisplayName', "Nutrient limitation")
hold on
plot(sigma_min, -param.z, '.r', 'LineWidth', 2,'MarkerSize',10, 'DisplayName', "Minimum")
hold on
plot(sigma_end, -param.z, '.k', 'LineWidth', 2,'MarkerSize',10, 'DisplayName', "\sigma")
legend('Location','north')
ylabel('Depth [m]')
grid on; grid minor

%% -------------------------------------------- Steady state solutions

figure('Name',"Steady-state-solutions")
subplot(3,3,[1 3])
plot(P(end,:),-param.z,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',3,'DisplayName',"Phytoplankton")
grid on; grid minor


subplot(3,3,[4 6])
plot(N(end,:),-param.z,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName',"Nutrients")
grid on; grid minor
ylabel('Depth [m]')

subplot(3,3,[7 9])
plot(D(end,:),-param.z,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',3,'DisplayName',"Detritus")
grid on; grid minor
xlabel('mmol N m^{-3}')

%% ------------------------------------------------------- NPD Dynamics

figure('Name',"NPZD Dynamics")
subplot(3,3,[1,3])
surface(t,-param.z',P')
shading interp
%ylabel('Depth')
%xlabel('Hours')
%title('[mmol N/m^{3}]')
%clim([0 0.2])
grid on; grid minor
colorbar
% Remove x-axis tick labels and tick marks
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);

subplot(3,3,[4,6])
surface(t,-param.z',N')
shading interp
ylabel('Depth [m]')
%xlabel('Hours')
%title('Nutrients [mmol N/m^{3}]')
grid on; grid minor
colorbar
% Remove x-axis tick labels and tick marks
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);

subplot(3,3,[7,9])
surface(t,-param.z',D')
shading interp
%ylabel('Depth')
xlabel('Hours')
%title('Detritus [mmol N/m^{3}]')
clim([0 1])
grid on; grid minor
colorbar

% Link the x-axis of all subplots
linkaxes(findall(gcf, 'type', 'axes'), 'x');


%% ----------------------------- Grid sensitivity analysis 
figure('Name',"Grid Sensitivity Analysis");

% Define different values of param.dz
deltaz_values = [param.dz/2,param.dz,param.dz*2,param.dz*8,param.dz*15,param.dz*30];

for j = 1:length(deltaz_values)
    % Update param.dz
    param.dz = deltaz_values(j);

    % Recalculate the number of grid cells and cell interval
    param.n = floor(depth / param.dz);
    param.z = 0.5 * param.dz : param.dz : depth - 0.5 * param.dz;

    % Initialize arrays to store results
    P2 = zeros(length(t), param.n);
    N2 = zeros(length(t), param.n);
    D2 = zeros(length(t), param.n);
    P2(1,:) = 0.1;
    N2(:,1) = param.Nb/10;

    % Solve ODEs
    [t2, y2] = ode45(@deriv, [0:time_duration], [P2(1,:) N2(1,:) D2(1,:)],[], param);

    % Extract P, N, D for the current grid spacing
    P2 = y2(:, 1:param.n);
    N2 = y2(:, param.n + 1:2 * param.n);
    D2 = y2(:, 2 * param.n + 1:3 * param.n);

    % Plot the results with varying color intensity
    color_index = (j - 1) / (length(deltaz_values) - 1); % Gradient index
    color = [1-color_index, color_index, 0]; % RGB color based on the index
    plot(P2(end, :), -param.z, 'LineWidth', 2, 'Color', color, 'DisplayName', ['\Delta z = ', num2str(param.dz)]);
    drawnow
    hold on;
end

title('Grid Sensitivity Analysis');
xlabel('Phytoplankton Concentration [mmol N/m^3]');
ylabel('Depth [m]');
legend('Location', 'best');
grid on; grid minor;




% ------------------------ Parameter sensitivity analysis

depth = 300; % Depth of water column [m]
param.dz = 1; % Length or depth of each cell
param.n = depth / param.dz; % Number of grid cells
param.z = 0.5 * param.dz:param.dz:depth - 0.5 * param.dz; % Cell interval

% Define different values of gmax
gmax_values = [param.gmax, param.gmax/10, param.gmax*10, param.gmax*100];

% Define colormap
colormap('winter')

figure('Name', 'Sensitivity analysis')

% Iterate over each value of gmax
for i = 1:length(gmax_values)
    % Update gmax value in param struct
    param.gmax = gmax_values(i);

    % Ensure P, N, and D are properly sized
    P1 = zeros(1, param.n);
    N1 = zeros(1, param.n);
    D1 = zeros(1, param.n);

    % Set initial conditions
    initial_conditions = zeros(3 * param.n, 1); % Create a column vector
    initial_conditions(1:param.n) = P1; % Assign initial values of P
    P1(1, 1) = P(1, 1);
    initial_conditions(param.n + 1:2 * param.n) = N1; % Assign initial values of N
    initial_conditions(2 * param.n + 1:3 * param.n) = D1; % Assign initial values of D
    D1(1, 1) = D(1, 1);

    % Call ode23 solver with the correct initial conditions
    [t1, y1] = ode45(@deriv, [0:time_duration], [P1 N1 D1], [], param);

    % Extract P, N, D for the current value of gmax
    P1 = y1(:, 1:param.n);
    N1 = y1(:, param.n + 1:2 * param.n);
    D1 = y1(:, 2 * param.n + 1:3 * param.n);

    % Define color based on gradient index
    color_index = (i - 1) / (length(gmax_values) - 1); % Gradient index
    color = [1-color_index, color_index, 0]; % RGB color based on the index

    % Plot P at the end of the simulation with gradient color
    plot(P1(end, :), -param.z, 'MarkerSize', 15, 'LineWidth', 2.5, 'Color', color, 'DisplayName', ["g_{max} = " + num2str(gmax_values(i)) + " h^{-1}"]);
    drawnow
    hold on;
end
legend('Location', 'southeast')
grid on; grid minor
title('[mmol N/m^{3}]')
xlabel('Phytoplankton [mmol N/m^{3}]')
ylabel('Depth [m]')



function I=calclight(t,P,D,param)
light=param.I0-param.I0.*sin(2*pi/(24*365).*t);
integral=cumsum((param.k_w+param.k_c.*(P+D)).*param.dz+1/2*param.dz.*(param.k_w+param.k_c.*(P+D)));
%I=light'.*exp(-param.k_w.*param.z-param.k_c.*integral'); %I=param.I0*exp(-integral');
I=light'.*exp(-integral');
end

function dydt = deriv(t,y,param) %calculates the dynamics of the inputs throughout the water column

P=y(1:param.n);
N=y(param.n+1:2*param.n);
D=y(2*param.n+1:3*param.n);

ix=2:param.n;

JaP(ix)=param.u*P(ix-1);
JdP(ix)=-param.D*(P(ix)-P(ix-1))/param.dz;

JaN(ix)=0;
JdN(ix)=-param.D*(N(ix)-N(ix-1))/param.dz;

JaD(ix)=-param.w*(D(ix-1)); 
JdD(ix)=-param.D*(D(ix)-D(ix-1))/param.dz;


%Phytoplankton
JaP(1)=0;
JaP(param.n+1)=0;
JdP(1)=0;
JdP(param.n+1)=0;

%Nutrients
JaN(1)=0;
JaN(param.n+1)=0; 
JdN(1)=0;
JdN(param.n+1)=-param.D*(param.Nb-N(end))/param.dz; 

%Detritus
JaD(1)=0;
JaD(param.n+1)=0; 
JdD(1)=0;
JdD(param.n+1)=-param.w*D(param.n); %0 istedet for måske


JP=JaP+JdP;
JN=JaN+JdN;
JD=JaD+JdD;

%Assembling the equations
I = calclight(t,P,D,param);

sigma_N=(N./(param.k_N+N));
sigma_L=(param.alpha.*I)./(sqrt(param.gmax^2+param.alpha^2.*I.^2));
sigma=(sigma_L'.*sigma_N);
%sigma_min=min(sigma_L',sigma_N);


gradientP=(-(JaP(2:param.n+1)+JdP(2:param.n+1)-JaP(1:param.n)-JdP(1:param.n))/param.dz)';
gradientN=(-(JaN(2:param.n+1)+JdN(2:param.n+1)-JaN(1:param.n)-JdN(1:param.n))/param.dz)';
gradientD=(-(JaD(2:param.n+1)+JdD(2:param.n+1)-JaD(1:param.n)-JdD(1:param.n))/param.dz)';


dPdt = param.gmax .* sigma .* P - param.eps .* P - param.m .* (P .* P) + gradientP; %new equation utilizing equations from the paper
dNdt = -param.gmax .* sigma .* P + param.thao .* D + gradientN; %new equation utilizing equations from the paper
dDdt = param.eps .* P + param.m .* (P .* P) - param.thao .* D + gradientD; %new equation utilizing equations from the paper


dydt=[dPdt;dNdt;dDdt];

end