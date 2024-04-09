function [dYdt,pI,pN] = odefun4(t,Y,param)
% Setup parameters
n = param.n; % no. cells
dz = param.dz;
%z = param.z;
u = param.u;
d = param.d;
l = param.l;
HI = param.HI;
HN = param.HN;
pmax = param.pmax;
Iin = param.Iin;
k = param.k;
Kbg = param.Kbg;
y = param.y;
ND = param.ND;
gamma = param.gamma;
r = param.r;
w = param.w;

% Splitting the Y variable into P (plankton) and N (nutrients):
P = Y(1:n);
N = Y((n+1):(2*n));
D = Y((2*n+1):end);

% Advection+diffusion flux for plankton:
Ja = zeros(n+1,1);
Jd = zeros(n+1,1);
for i = 2:n
    Ja(i) = u*P(i-1);
    Jd(i) = -d*(P(i)-P(i-1))/dz;
    
end
JP = Ja+Jd; % Total plankton flux
% Boundary conditons:
JP(1) = 0;
JP(end) = 0;

% Diffusion flux for nutrients:
JN = zeros(n+1,1);
for i = 2:n
    JN(i) = -d*(N(i)-N(i-1))/dz;
    
end
% Boundary conditions:
JN(1) = 0;
JN(end) = -d*(ND-N(n-1))/dz;


% Advection+diffusion flux for detritus:
JaD = zeros(n+1,1);
JdD = zeros(n+1,1);
for i = 2:n
    JaD(i) = w*D(i-1);
    JdD(i) = -d*(D(i)-D(i-1))/dz;
    
end
JD = JaD+JdD; % Total plankton flux
% Boundary conditons:
JD(1) = 0;
JD(end) = w*D(end);

% Light intensity as a function of depth:
I = Iin*exp(-cumsum( (P*k + Kbg)*dz));

for i = 1:n
    
    % Production due to light or minerals (limiting factor):
    pI(i) = I(i)/(HI+I(i));
    pN(i) = N(i)/(HN+N(i));
    p(i) = pmax*min( pI(i) , pN(i) ); % growth due to minerals or light
    %g(i) = p(i) - l; % growth minus losses
    
    % PDEs
    dPdt(i) = p(i)*P(i) - l*P(i) - gamma*P(i)^2*y - (JP(i+1)-JP(i))/dz;      % Plankton [cells/m3]
    %dPdt(i) = p(i)*P(i) - l*P(i) - (JP(i+1)-JP(i))/dz;
    dNdt(i) = - y*p(i)*P(i) + r*D(i) - (JN(i+1)-JN(i))/dz;                 % Nutrient [mmol nutrient/m3]
    dDdt(i) = l*P(i)*y + gamma*(P(i)*y)^2 - r*D(i) - (JD(i+1)-JD(i))/dz;       % Detritus [mmol detritus/m3]
    %dDdt(i) = l*P(i)*y - r*D(i) - (JD(i+1)-JD(i))/dz;
    
end
pI = pI';
pN = pN';

dPdt = dPdt';
dNdt = dNdt';
dDdt = dDdt';
dYdt = [dPdt;dNdt;dDdt];

end