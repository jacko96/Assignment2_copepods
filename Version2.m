clear all;clc;
close all

%System parameters
depth=300; %depth of water column [m]
param.dz=1; %length or depth of each cell
param.n=depth/param.dz; %number of gridcells
param.z=0.5*param.dz:param.dz:depth-0.5*param.dz; %cell interval
time_duration=1000; %[hours] 

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
param.m=1.5/24*14; % grazing [m^3 (mg N)^-1 h^-1]
param.thao=0.1/24;%0.01; % remineralization[h^-1]
%param.Av=5*10^(-5)*(60*60); % %Vertical diffusion [m^2 h^-1]
param.w=15/24; % sinking of detritus [m h^-1]
param.u=0.001/24; %vertical velocity of phytoplankton [m h^-1]
param.Nb=30*14; %bottom concentration of nutrients [mg N m^-3]


%Copepods Parameters
param.S=5; %Number of lifestages for Copepods

param.CtoN=5.6;
param.v=0.0052/24*param.CtoN*0.001; %passive clearance rate coefficient [L mg N^{-3/4} h^{-1}]
param.q=-1/4; %passive clearance rate exponent (no exponent)
param.h=0.4*param.CtoN/24*0.001; %maximum ingestion rate coefficent [mg N^{1/4} h^{-1}]
param.h_n=-1/4; %maximum ingestion rate exponent
param.kappa=0.048/24*param.CtoN*0.001; % respiration rate coefficent [mg N^{1/4} h^{-1}]
param.p=-1/4; %respiration rate exponent
param.assi=0.67; %assimilation rate
param.reproc=0.25; %reproduction efficiency rate
param.mu_ht=0.003*param.CtoN*0.001; %Mortality by higher trophic levels [mg N^{1/4} mgN^{-2} L^{-2} h^{-1}]
param.Cs=0.001; %Copepod speed velocity [m/h]
param.b=0.37/24; %reproduction rate [h^{-1}]
param.mu=0.01; %death or mortality rate of copepods [h^-1] --> based on averages found online

%Copepods grid space
param.m0=0.01; % mgN pr individual
param.ma=1; % mgN pr individual
param.m_minus=exp(linspace(log(param.m0),log(param.ma),param.S+1)); % mgN pr individual
param.m_plus=exp(((log(param.m_minus(2:param.S+1))+log(param.m_minus(1:param.S)))/2)); % mgN pr individual
param.m_ratio=param.m_plus(2)./param.m_plus(1);


%--------------------------------- Initialization


%C1=0.0*param.z;
C1=ones(1,param.n);
%C2=0.0*param.z;
%C2(1)=1;
C2=ones(1,param.n);
%C3=0.0*param.z;
%C3(1)=1;
C3=ones(1,param.n);
%C4=0.0*param.z;
%C4(1)=1;
C4=ones(1,param.n);
%C5=0.0*param.z;
%C5(1)=1;
C5=ones(1,param.n);

%C1=C1*0;
%C2=C2*0;
%C3=C3*0;
%C4=C4*0;
%C5=C5*0;

%P = 1 * param.z; %
P = 50*ones(1,param.n);
%N = 1 * param.z; %
N = 100*ones(1,param.n);
D = 0 * param.z; %

N(1:end)=0+param.Nb/length(param.z)*param.z;
P(1,1)=0.1;

Initial=[P N D C1 C2 C3 C4 C5];


%options1 = odeset('Refine',1);
%options2 = odeset(options1,'NonNegative',1);
[t,y]=ode45(@deriv,[0:time_duration],Initial,[],param);


P = y(:,1:param.n);
N = y(:,param.n+1:2*param.n);
D = y(:,2*param.n+1:3*param.n);

C1 =y(:,3*param.n+1:4*param.n);
C2 =y(:,4*param.n+1:5*param.n);
C3 =y(:,5*param.n+1:6*param.n);
C4 =y(:,6*param.n+1:7*param.n);
C5 =y(:,7*param.n+1:8*param.n);


%% ------------------------------------------------------- NPD Dynamics

figure('Name', "NPZD Dynamics");

% Subplot for Phytoplankton (P)
subplot(3, 3, 1);
surface(t, -param.z', P');
shading interp;
title('Phytoplankton');
colorbar;
grid on; grid minor;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);

% Subplot for Nutrients (N)
subplot(3, 3, 2);
surface(t, -param.z', N');
shading interp;
title('Nutrients');
colorbar;
grid on; grid minor;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);

% Subplot for Detritus (D)
subplot(3, 3, 3);
surface(t, -param.z', D');
shading interp;
title('Detritus');
colorbar;
grid on; grid minor;
set(gca, 'XTickLabel', []);
set(gca, 'XTick', []);

% Subplot for C1
subplot(3, 3, 4);
surface(t, -param.z', C1');
shading interp;
title('C1');
colorbar;
grid on; grid minor;
ylabel('Depth [m]');

% Subplot for C2
subplot(3, 3, 5);
surface(t, -param.z', C2');
shading interp;
title('C2');
colorbar;
grid on; grid minor;

% Subplot for C3
subplot(3, 3, 6);
surface(t, -param.z', C3');
shading interp;
title('C3');
colorbar;
grid on; grid minor;
xlabel('Hours');

% Subplot for C4
subplot(3, 3, 7);
surface(t, -param.z', C4');
shading interp;
title('C4');
colorbar;
grid on; grid minor;
xlabel('Hours');
ylabel('Depth [m]');

% Subplot for C5
subplot(3, 3, 8);
surface(t, -param.z', C5');
shading interp;
title('C5');
colorbar;
grid on; grid minor;
xlabel('Hours');

% Adjust the layout
sgtitle('State Variables Dynamics [mgN]');
linkaxes(findall(gcf, 'type', 'axes'), 'x');


%%
figure()
surface(t, -param.z', C1');
shading interp;
title('C2');
colorbar;
grid on; grid minor;

%% ----------------------------- Grid sensitivity analysis 
% figure('Name',"Grid Sensitivity Analysis");
% 
% % Define different values of param.dz
% deltaz_values = [param.dz/2,param.dz,param.dz*2,param.dz*8,param.dz*15,param.dz*30];
% 
% for j = 1:length(deltaz_values)
%     % Update param.dz
%     param.dz = deltaz_values(j);
% 
%     % Recalculate the number of grid cells and cell interval
%     param.n = floor(depth / param.dz);
%     param.z = 0.5 * param.dz : param.dz : depth - 0.5 * param.dz;
% 
%     % Initialize arrays to store results
%     P2 = zeros(length(t), param.n);
%     N2 = zeros(length(t), param.n);
%     D2 = zeros(length(t), param.n);
%     P2(1,:) = 0.1;
%     N2(:,1) = param.Nb/10;
% 
%     % Solve ODEs
%     [t2, y2] = ode45(@deriv, [0:time_duration], [P2(1,:) N2(1,:) D2(1,:)],[], param);
% 
%     % Extract P, N, D for the current grid spacing
%     P2 = y2(:, 1:param.n);
%     N2 = y2(:, param.n + 1:2 * param.n);
%     D2 = y2(:, 2 * param.n + 1:3 * param.n);
% 
%     % Plot the results with varying color intensity
%     color_index = (j - 1) / (length(deltaz_values) - 1); % Gradient index
%     color = [1-color_index, color_index, 0]; % RGB color based on the index
%     plot(P2(end, :), -param.z, 'LineWidth', 2, 'Color', color, 'DisplayName', ['\Delta z = ', num2str(param.dz)]);
%     drawnow
%     hold on;
% end
% 
% title('Grid Sensitivity Analysis');
% xlabel('Phytoplankton Concentration [mmol N/m^3]');
% ylabel('Depth [m]');
% legend('Location', 'best');
% grid on; grid minor;




%% ------------------------ Parameter sensitivity analysis

% depth = 300; % Depth of water column [m]
% param.dz = 1; % Length or depth of each cell
% param.n = depth / param.dz; % Number of grid cells
% param.z = 0.5 * param.dz:param.dz:depth - 0.5 * param.dz; % Cell interval
% 
% % Define different values of gmax
% gmax_values = [param.gmax, param.gmax/10, param.gmax*10, param.gmax*100];
% 
% % Define colormap
% colormap('winter')
% 
% figure('Name', 'Sensitivity analysis')
% 
% % Iterate over each value of gmax
% for i = 1:length(gmax_values)
%     % Update gmax value in param struct
%     param.gmax = gmax_values(i);
% 
%     % Ensure P, N, and D are properly sized
%     P1 = zeros(1, param.n);
%     N1 = zeros(1, param.n);
%     D1 = zeros(1, param.n);
% 
%     % Set initial conditions
%     initial_conditions = zeros(3 * param.n, 1); % Create a column vector
%     initial_conditions(1:param.n) = P1; % Assign initial values of P
%     P1(1, 1) = P(1, 1);
%     initial_conditions(param.n + 1:2 * param.n) = N1; % Assign initial values of N
%     initial_conditions(2 * param.n + 1:3 * param.n) = D1; % Assign initial values of D
%     D1(1, 1) = D(1, 1);
% 
%     % Call ode23 solver with the correct initial conditions
%     [t1, y1] = ode45(@deriv, [0:time_duration], [P1 N1 D1], [], param);
% 
%     % Extract P, N, D for the current value of gmax
%     P1 = y1(:, 1:param.n);
%     N1 = y1(:, param.n + 1:2 * param.n);
%     D1 = y1(:, 2 * param.n + 1:3 * param.n);
% 
%     % Define color based on gradient index
%     color_index = (i - 1) / (length(gmax_values) - 1); % Gradient index
%     color = [1-color_index, color_index, 0]; % RGB color based on the index
% 
%     % Plot P at the end of the simulation with gradient color
%     plot(P1(end, :), -param.z, 'MarkerSize', 15, 'LineWidth', 2.5, 'Color', color, 'DisplayName', ["g_{max} = " + num2str(gmax_values(i)) + " h^{-1}"]);
%     drawnow
%     hold on;
% end
% legend('Location', 'southeast')
% grid on; grid minor
% title('[mmol N/m^{3}]')
% xlabel('Phytoplankton [mmol N/m^{3}]')
% ylabel('Depth [m]')



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


% ---------------------------- Fluxes

ix=2:param.n; %Vector notation

%Phytoplankton
JaP(ix)=param.u*P(ix-1);
JdP(ix)=-param.D*(P(ix)-P(ix-1))/param.dz;

%Nutrients
JaN(ix)=0;
JdN(ix)=-param.D*(N(ix)-N(ix-1))/param.dz;

%Detritus
JaD(ix)=-param.w*(D(ix-1)); 
JdD(ix)=-param.D*(D(ix)-D(ix-1))/param.dz;

%Copepods
JaC(ix)=0;
JdC(ix)=-param.Cs*(P(ix)-P(ix-1))/param.dz; %the movement of copepods toward the phytoplankton

% --------------------------------- Boundary conditions
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

%Copepods
JaC(1)=0;
JaC(param.n+1)=0;
JdC(1)=0;
JdC(param.n+1)=0;


JP=JaP+JdP;
JN=JaN+JdN;
JD=JaD+JdD;
JC=JaC+JdC;

%Assembling the equations
I = calclight(t,P,D,param);

sigma_N=(N./(param.k_N+N));
sigma_L=(param.alpha.*I)./(sqrt(param.gmax^2+param.alpha^2.*I.^2));
sigma=(sigma_L'.*sigma_N);
%sigma_min=min(sigma_L',sigma_N);

% ---------------------------- Advection-Diffusion
gradientP=(-(JaP(2:param.n+1)+JdP(2:param.n+1)-JaP(1:param.n)-JdP(1:param.n))/param.dz)';
gradientN=(-(JaN(2:param.n+1)+JdN(2:param.n+1)-JaN(1:param.n)-JdN(1:param.n))/param.dz)';
gradientD=(-(JaD(2:param.n+1)+JdD(2:param.n+1)-JaD(1:param.n)-JdD(1:param.n))/param.dz)';
gradientC=(-(JaC(2:param.n+1)+JdC(2:param.n+1)-JaC(1:param.n)-JdC(1:param.n))/param.dz)';

% ----------------------------- Copepods

% Feeding level (0-1):
%f=(param.v.*(param.m_plus.^param.q).*P)./(param.v.*(param.m_plus.^param.q).*P + param.h.*param.m_plus.^param.h_n);
f = (param.v.*(param.m_plus.^param.q).*P/1000)./(param.v.*(param.m_plus.^param.q).*P/1000 + param.h.*param.m_plus.^param.h_n); % P/1000 fordi der er 1000 L på en m3.

fc = param.kappa./(param.assi*param.h).*param.m_plus.^(param.p - param.h_n);
% biomass accumulation/production rate:
v = param.assi*param.h.*param.m_plus.^param.h_n.*(f - fc);
g = max(0,v); % net energy gain
mu_st = min(0,v); %mu constant?
mu = mu_st + param.mu; % copepod mortality
b = param.reproc*g(end,end); % copepod birth rate
gamma = (g - mu)./(1 - (param.m_ratio).^(1 - mu./g)); % Overgang fra lille til stor copepod
gamma = max(zeros(size(gamma)),gamma); % hvis gamma bliver negativ

dC1dt = b.*C5 + g(1)*C1 - gamma(1)*C1 - mu(1)*C1 + gradientC;
dC2dt = gamma(1)*C1 + g(2)*C2 - gamma(2)*C2 - mu(2)*C2 + gradientC;
dC3dt = gamma(2)*C2 + g(3)*C3 - gamma(3)*C3 - mu(3)*C3 + gradientC;
dC4dt = gamma(3)*C3 + g(4)*C4 - gamma(4)*C4 - mu(4)*C4 + gradientC;
dC5dt = gamma(4)*C4 - mu(5)*C5 + gradientC;

%dC1dt=gradientC;
%dC2dt=gradientC;
%dC3dt=gradientC;
%dC4dt=gradientC;
%dC5dt=gradientC;

% ----------------------------- NPD equations
%dPdt = param.gmax .* sigma .* P - param.eps .* P - param.m .* (C1+C2+C3+C4+C5) + gradientP; %new equation utilizing equations from the paper
dPdt = param.gmax .* sigma .* P - param.eps .* P - sum(param.v.*(param.m_plus.^param.q).*P/1000,2) + gradientP; % Har erstattet grazing rate med encounter rate fra ligning (1)
dNdt = -param.gmax .* sigma .* P + param.thao .* D + gradientN; %new equation utilizing equations from the paper
dDdt = param.eps .* P + param.m .* (P .* P) - param.thao .* D +param.thao*(mu(1)*C1+mu(2)*C2+mu(3)*C3+mu(4)*C4+mu(5)*C5)+gradientD; %new equation utilizing equations from the paper

%adding the loss of copepods to detritus equation



dydt=[dPdt;dNdt;dDdt;dC1dt;dC2dt;dC3dt;dC4dt;dC5dt];

end