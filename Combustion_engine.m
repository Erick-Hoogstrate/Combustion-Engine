clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. Ready for use
global Runiv Pref
Runiv=8.314472;
Pamb=1.01235e5; % Reference pressure, 1 atm!
Tamb=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. For the final assignment take the ones from the specific case you are supposed to do.                           
cFuel='Gasoline';   
%% Select all species
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important these are molefractions
MAir = Xair*Mi';
Yair = Xair.*Mi/MAir;
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T

%% ---------------Code------------
%% Stoichiometric combustion
% (100-percentage)*Cx_gasHy_gasOz_gas + percentage*Cx_ethHy_ethOz_eth +
% ((x_gas+x_eth)+(y_gas+y_eth)/4-(z_gas+z_eth)/2)(O2 + 3.76N2) -->
% (x_gas+x_eth)CO2 + ((y_gas+y_eth)/2)H2O + ((x_gas+x_eth)+(y_gas+y_eth)/4-(z_gas+z_eth)/2)N2

percentage = 5;

% mol_gas = (100-percentage)*volume*density * molaire massa
mol_gas = (100-percentage)*1*0.75* SpS(1).Mass*1000;                     % Density gasoline is from 0.71 to 0.77 g/cm3

% mol_ethanol = percentage*volume*density * molaire massa
mol_ethanol = percentage*1*0.78945* 46;                                  % Density Ethanol    0.78945 g/cm3 (at 20 ?C)


x_gas = mol_gas*SpS(1).Elcomp(3); 
y_gas = mol_gas*SpS(1).Elcomp(2); 
z_gas = mol_gas*SpS(1).Elcomp(1); 

x_eth = mol_ethanol*2;
y_eth = mol_ethanol*6;
z_eth = mol_ethanol*1;

%% Molair masses
C = Sp(myfind({Sp.Name},{'C'})).Mass * 1000;      % Carbon    [g/mol]
H = Sp(myfind({Sp.Name},{'H'})).Mass * 1000;      % Hydrogen  [g/mol]
O = Sp(myfind({Sp.Name},{'O'})).Mass * 1000;      % Oxygen    [g/mol]
N = Sp(myfind({Sp.Name},{'N'})).Mass * 1000;      % Nitrogen  [g/mol]
%% AF Gasoline blend
M_gas = C*x_gas + H*y_gas + O*z_gas;                        % [g]
M_eth = C*x_eth + H*y_eth + O*z_eth;                        % [g]

M_air = ((x_gas+x_eth)+(y_gas+y_eth)/4-(z_gas+z_eth)/2)*(O*2 + 3.76*N*2);   % [g]
% AF_stoi = (x_gas+y_gas/4-z_gas/2)*(O*2 + 3.76*N*2)/M_gas
AF_stoi = M_air/(M_gas+M_eth)
% lambda = AF/AF_stoi
%%
%state variables
p1 = Pamb;
T1 = Tamb;

%engine geometric parameters
bore = 0.068; %m
stroke = 0.054; %m
radius = bore/2; %m
rod = 0.091313; %m
r = 8.5;                                                                    %compression ratio
RPM = 3000;

%calculating swept volume and the clearance volume
v_d = (pi/4)*bore^2*stroke; %m^3
v_c = v_d/(r-1); %m^3


%point 1
v1 = v_d + v_c; %m^3
v2 = v_c; %m^3
g = 1.4; %gamma

%starting point / point 0
p0 = p1;
v0 = (v_d + v_c)/r
T0 =(p0*v0*T1/(p1*v1));

%calculation of state variables at state point 2
%p2*v2^g=p1*v1^g
p2 = p1*r^g;
c = p1*v1^g;
v_compression = engine_kinematics(bore,stroke,r,rod,180,0);
p_compression =(c./v_compression.^g);
%(p1*v1/T1)=(p2*v2/T2)|T2=(p2*v2*T2)/(p1*v1) 
T2 =(p2*v2*T1/(p1*v1));

%calculation of state varibles at point 3 
%(p2/T2)=(p3/T3)|v3=v2
%input
hv_gas = 45*10^6 %J/kg. Between 44-46. This is de heating value for gasoline (value from the internet)
hv_eth = 29.7*10^6 %J/kg. This is de heating value for ethanol
T3 = 0.5*(AF_stoi*10^-3)*(((hv_gas*((M_gas*10^-3)))+(hv_eth*(M_eth*10^-3)))/((cv_function(percentage,T2))*((M_gas*10^-3)+(M_eth*10^-3))))%Calculated with heating value [J/kg], the stoichiometric air to fuel ratio (AF_stoi), the specific heat capacity of fuel and air [J/kg K], the air and fuel inlet temperatures [K]
%T3=998; %between 923-1073 Kelvin
p3 = (p2*T3/T2);
v3 = v2;

%calculation of state varibes at point 4
%p3*v3^gamma=p4*v4^gammma
v4 = v1;
p4 = p3*(v3/v4)^g;
c2 = p3*(v3)^g;
v_expansion = engine_kinematics(bore,stroke,r,rod,0,180);
p_expansion = c2./(v_expansion.^g);


%plot figure
figure(1)
hold on
plot(v0,p0,'*','color','r')
plot(v1,p1,'*','color','r')
plot(v2,p2,'*','color','r')
plot(v3,p3,'*','color','r')
plot(v4,p4,'*','color','r')
set(gca,'FontSize',30)
P0=plot([v0 v1],[p0 p1],'-','color','blue','LineWidth',8);%isobaric expansion
P1=plot(v_compression,p_compression,'color','m','LineWidth',8);%isentropic compression
P2=plot([v2 v3],[p2 p3],'color','b','LineWidth',8);%isochoric heat addition
P3=plot(v_expansion,p_expansion,'color','k','LineWidth',8);%isentropic expansion
P4=plot([v4 v1],[p4,p1],'color','g','LineWidth',8);%isochoric heat rejection
P5=plot([v1 v0],[p1 p0],'--','color','yellow','LineWidth',8);%isobaric compresion
xlabel('Volume [m^3]')
ylabel('Pressure [Pa]')
title({['PV Diagram'],[(sprintf('Otto cycle E%.0f', percentage))]}, 'FontSize', 50)
legend([P0 P1 P2 P3 P4 P5],{'isobaric expansion','isentropic compression','isochoric heat addition','isentropic expansion','isochoric heat rejection', 'isobaric compression'})

%thermal efficiency
thermalefficiency=(1-1/r^(g-1))    

%%
%Work
area12 = trapz(-v_compression,p_compression)                                 %start for calculation work
area34 = trapz(v_expansion,p_expansion)

work = area34-area12 %pa*m^3 is equal to J. Work per cycle

% %%
% %Heat losses
% 
% %constants
% Tw = 323;      % Estimation of the wall temperature of the engine (50 degrees celcius)
% C0 = 120;    % between 110-130
% % P =       % dependent on crank angle function
% V_mp = 2*stroke*(RPM/60);
% P_mot = 5*Pamb;
% 
% % h_g = C0*bore^-0.2*P^0.8*((C1*V_mp)+(C2*(V_d*T1)/(p1*v1)*(P-P_mot)))^0.8*T^-0.53
% 
% % step 0-1 (isobaric expansion)
% C1 = 6.18;
% C2 = 0;
% P_begin = p0;
% P_end = p1;
% T_begin = T0;
% T_end = T1;
% V_begin = v0;
% V_end = v1;
% % A_01 = 2*v1/radius  
% 
% for T=T_begin:T_end
%     for P=P_begin:P_end
%         h_g_01 = C0*bore^-0.2*P^0.8*((C1*V_mp))^0.8*T^-0.53;
%     end
%     for v=V_begin:V_end
%         A_01 = 2*v/radius;
%     end
%     dQdt = -h_g_01*A_01*(T-Tw);
% end
% 
% 
% % step 1-2 (isentropic compression)
% C1 = 2.28;
% C2 = 0;
% P = p_compression;
% % P_begin = p0;
% % P_end = p1;
% T_begin = T1;
% T_end = T2;
% A_12 = 2*v_compression/radius;  
% % T = (p_compression.*v_compression*T1/(p1*v1);     % Unclear what T should be used and whether this should be a constant value or a function
% 
% for T=T_begin:T_end
% %     for P=P_begin:P_end
%     h_g_12 = C0*bore^-0.2*P.^0.8.*((C1*V_mp))^0.8*T^-0.53;
% %     end
%     dQdt = -h_g_12.*A_12.*(T-Tw);
% end
% 
% 
% % step 2-3 (isochoric heat addition)
% C1 = 2.28;
% C2 = 3.24*10^-3;
% P_begin = p2;
% P_end = p3;
% T_begin = T2;
% T_end = T3;
% A_23 = 2*v2/radius;  
% 
% for T=T_begin:T_end
%     for P=P_begin:P_end
%         h_g_23 = C0*bore^-0.2*P^0.8*((C1*V_mp)+(C2*(v_d*T1)/(p1*v1)*(P-P_mot)))^0.8*T^-0.53;
%     end
%     dQdt = -h_g_23*A_23*(T-Tw);
% end
% 
% 
% % step 3-4 (isentropic expansion)
% C1 = 2.28;
% C2 = 3.24*10^-3;
% P = p_expansion;
% % P_begin = p3;
% % P_end = p4;
% T_begin = T3;
% T_end = T4;
% A_34 = 2*v_expansion/radius;  
% 
% for T=T_begin:T_end
% %     for P=P_begin:P_end
%     h_g_34 = C0*bore^-0.2*P^0.8*((C1*V_mp)+(C2*(v_d*T1)/(p1*v1)*(P-P_mot)))^0.8*T^-0.53;
% %     end
%     dQdt = -h_g_01*A_34*(T-Tw);
% end
% 
% 
% % Unknown which coefficients need to be used
% % % step 4-1 (isochoric heat rejection)
% % C1 = 6.18;
% % C2 = 0;
% % P = p0;
% % T = T1;
% % A_41 = 2*v1/radius;  
% % 
% % h_g_01 = C0*bore^-0.2*P^0.8*((C1*V_mp))^0.8*T^-0.53;
% % dQdt = -h_g_01*A_41*(T-Tw);
% 
% 
% % step 1-0 (isobaric compression)
% C1 = 6.18;
% C2 = 0;
% P_begin = p1;
% P_end = p0;
% T_begin = T1;
% T_end = T0;
% V_begin = v1;
% V_end = v0;
% % A_10 = 2*v/radius;  
% 
% for T=T_begin:T_end
%     for P=P_begin:P_end
%         h_g_10 = C0*bore^-0.2*P^0.8*((C1*V_mp))^0.8*T^-0.53;
%     end
%     for v=V_begin:V_end
%         A_10 = 2*v/radius; 
%     end
%     dQdt = -h_g_10*A_10*(T-Tw);
% end
% 
% 
% 
% return
% cv = cv_function(percentage,T)%cv based on percentage ethanol and temperature
% T4=(p4/p1)*T1%Calculated with ideal gas law
% Qc = cv*(MAir+MFuel)*(T4-T1)
% 
% %efficiency
% efficiency = work/(work+Qc)
% 
% 
% 
% return
%%
%Power
dt_no = 0.0378%Time per cycle, no load
dt_full = 0.040065%Time per cycle, full load
dt_half = 0.0389344545%Time per cycle, half load

power = work/dt_no%J/s is equal to W. Power per cycle

return

%%
% Initialisation
p(1)=P0;T(1)=T0;
pad(1)=p(1);
v_c = v_d/(r-1);
Ca(1)=61;
V(1)=v_c %(Ca(1),stroke,bore,rod,radius); % Vcyl is a function that computes cyl vol as fie of crank-angle for given B,S,l and rc
%m(1)=p(1)*V(1)/Runiv/T(1); %ye

m_E0_NL= 0.0000079; % Mass is constant, valves are closed, density=0.74 g/cm^3, mass per cycle
mf=m_E0_NL
m(1) = m_E0_NL;

a = 5
n = 3
%
%
% Loop over crank-angle, with 'for' construction
NCa=360; % Number of crank-angles
dCa=0.5; % Stepsize
NSteps=NCa/dCa;


for i=2:NSteps,
    Ca(i)=Ca(i-1)+dCa;
    V(i)= pi*(bore/2)^2*(rod+r-(r*cosd(Ca(i))+sqrt(rod^2-r^2*(sind(Ca(i)))^2)))+v_c; % New volume for current crank-angle
    m(i)=m(i-1); % Mass is constant, valves are closed
    dV=V(i)-V(i-1); % Volume change
    Q_LHV= 43.4e6
    theta_d= 240
    xb= 1-exp(-a*((Ca(i)-Ca(1)/theta_d)^n));
    dQcom = Q_LHV*mf*n*a*(1-xb)/theta_d*(Ca(i)-Ca(1)/theta_d)^(n-1); % Heat Release by combustion
    dT=(dQcom-p(i-1)*dV)/Cv/m(i-1); % 1st Law dU=dQ-pdV (closed system)
    A(i)= (pi/2)*bore^2+pi*bore*(stroke/2)*((1/r)+1-cosd(Ca(i))+sqrt((1/r)^2-sind(Ca(i)))^2);
    T(i)=T(i-1)+dT;
    p(i)=m(i)*Rg*T(i)/V(i); % Gaslaw
end;
