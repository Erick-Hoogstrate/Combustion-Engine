%%
% Initialisation
close all;

a = 5;
n = 3;
Q_LHV= 43.4e6;   % [J/kg]
theta_d= 60;     % ca difference start and end combustion 400-340
theta_s= 340;    % ca at start of combustion
Ca(1)=theta_s;

NCa = 360;       % Number of crank-angles
dCa = 0.5;       % Stepsize
NSteps = NCa/dCa;

m_E0_NL = 0.0000079; % Mass is constant, valves are closed, density=0.74 g/cm^3, mass per cycle
mf = m_E0_NL;
m(1) = m_E0_NL;


if Ca < theta_s;
   Ca - theta_s  == 0;
end

for i=2:NSteps,
    Ca(i)=Ca(i-1)+dCa;
    xb(i) = 1 - exp(-a*((Ca(i)-theta_s)/theta_d)^n);
    dQcom(i) = Q_LHV*mf*n*a*(1-xb(i))/theta_d*((Ca(i)-theta_s)/theta_d)^(n-1); % Heat Release by combustion
end;


figure()
plot(Ca,dQcom)
xlim([260 480])
xlabel('Crank angle (\theta)')
ylabel('Combustion heat release (J)')
title({'Wiebe function';'E0 No load'})