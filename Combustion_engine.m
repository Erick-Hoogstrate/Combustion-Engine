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
load = "no"; %no,half,full load

% mol_gas = (100-percentage)*volume*density * molaire massa
mol_gas = (100-percentage)*1*0.75/(SpS(1).Mass*1000);                     % Density gasoline is from 0.71 to 0.77 g/cm3

% mol_ethanol = percentage*volume*density * molaire massa
mol_ethanol = percentage*1*0.78945/46;                                  % Density Ethanol    0.78945 g/cm3 (at 20 ?C)


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
AF_stoi = M_air/(M_gas+M_eth);
% lambda = AF/AF_stoi
%%
%state variables
p1 = Pamb;
T1 = Tamb;

%engine geometric parameters
bore = 0.068; %m
stroke = 0.054; %m
radius = bore/2; %m
rod  = 0.091313;    %[m]
r_crank = stroke/2;       % Radius of the crankshaft[m] (half the stroke)
r_cyl   = bore/2;         % Radius cyllinder [m]
r = 8.5;                                                                    %compression ratio
RPM = 3000;

%calculating swept volume and the clearance volume
v_d = (pi/4)*bore^2*stroke; %m^3
v_c = v_d/(r-1); %m^3


%point 1
v1 = v_d + v_c; %m^3
v2 = v_c; %m^3
gamma = 1.4; %gamma

%starting point / point 0
p0 = p1;
v0 = (v_d + v_c)/r;
T0 =(p0*v0*T1/(p1*v1));

%calculation of state variables at state point 2
%p2*v2^g=p1*v1^g
p2 = p1*r^gamma;
c = p1*v1^gamma;
v_compression = engine_kinematics(bore,stroke,r,rod,180,0);
p_compression =(c./v_compression.^gamma);
%(p1*v1/T1)=(p2*v2/T2)|T2=(p2*v2*T2)/(p1*v1) 
T2 =(p2*v2*T1/(p1*v1));

%calculation of state varibles at point 3 
%(p2/T2)=(p3/T3)|v3=v2
%input
hv_gas = 45*10^6; %J/kg. Between 44-46. This is de heating value for gasoline (value from the internet)
hv_eth = 29.7*10^6; %J/kg. This is de heating value for ethanol
T3 = 0.5*(AF_stoi*10^-3)*(((hv_gas*((M_gas*10^-3)))+(hv_eth*(M_eth*10^-3)))/((cv_function(percentage,T2))*((M_gas*10^-3)+(M_eth*10^-3))));%Calculated with heating value [J/kg], the stoichiometric air to fuel ratio (AF_stoi), the specific heat capacity of fuel and air [J/kg K], the air and fuel inlet temperatures [K]
%T3=998; %between 923-1073 Kelvin
if load == "no"
    p3 = (p2*T3/T2)*0.1;
elseif load == "half"
    p3 = (p2*T3/T2)*0.15;
elseif load == "full"
    p3 = (p2*T3/T2)*0.3;
end
v3 = v2;

%calculation of state varibes at point 4
%p3*v3^gamma=p4*v4^gammma
v4 = v1;
p4 = p3*(v3/v4)^gamma;
c2 = p3*(v3)^gamma;
v_expansion = engine_kinematics(bore,stroke,r,rod,0,180);
p_expansion = c2./(v_expansion.^gamma);


%plot figure
% figure(1)
% hold on
% plot(v0,p0,'*','color','r')
% plot(v1,p1,'*','color','r')
% plot(v2,p2,'*','color','r')
% plot(v3,p3,'*','color','r')
% plot(v4,p4,'*','color','r')
% set(gca,'FontSize',30)
% P0=plot([v0 v1],[p0 p1],'-','color','blue','LineWidth',8);%isobaric expansion
% P1=plot(v_compression,p_compression,'color','m','LineWidth',8);%isentropic compression
% P2=plot([v2 v3],[p2 p3],'color','b','LineWidth',8);%isochoric heat addition
% P3=plot(v_expansion,p_expansion,'color','k','LineWidth',8);%isentropic expansion
% P4=plot([v4 v1],[p4,p1],'color','g','LineWidth',8);%isochoric heat rejection
% P5=plot([v1 v0],[p1 p0],'--','color','yellow','LineWidth',8);%isobaric compresion
% xlabel('Volume [m^3]')
% ylabel('Pressure [Pa]')
% title({['PV Diagram'],(sprintf('Otto cycle E%.0f %s load', [percentage, load]))}, 'FontSize', 50)
% legend([P0 P1 P2 P3 P4 P5],{'isobaric expansion','isentropic compression','isochoric heat addition','isentropic expansion','isochoric heat rejection', 'isobaric compression'})

%thermal efficiency
thermalefficiency=(1-1/r^(gamma-1));    

%% REAL CYCLE
%% Intake pressure

intake_pressure = [0.21 0.34 0.44; 0.28 0.39 0.46; 0.31 0.41 0.50; 0.23 0.35 0.56].*10^5; % rows represent fuels (E0 upper row), columns represent loads (NL first column)
                                                                                           % values derived from average pressure plots experimental data
fuel_mass_exp   = [4.071 5.3597 7.9903; 3.3331 4.9078 6.7648; 3.8139 4.8079 6.0022; 3.1186 4.2123 5.8703].*10^-6;   % idem 
                                                                                     
if percentage == 0
    E_row = 1;
elseif percentage == 5
    E_row = 2;
elseif percentage == 10
    E_row = 3;
elseif percentage == 15
    E_row = 4;
end

if load == "no"
    load_column = 1;
elseif load == "half"
    load_column = 2;
elseif load == "full"
    load_column = 3;
end

if ismember(percentage,[0, 5, 10, 15])
    p_int = intake_pressure(E_row,load_column);     % intake pressure
    mf = fuel_mass_exp(E_row,load_column);          % fuel mass per cycle
else
    Ep = percentage;
    p_int = intake_pressure(mod(Ep,5),load_column)+ (intake_pressure((5-mod(Ep,5)),load_column)-intake_pressure(mod(Ep,5),load_column))*(Ep-(Ep-mod(Ep,5)))/((Ep+(5-mod(Ep,5)))-(Ep-mod(Ep,5)));     % intake pressure
    mf = fuel_mass_exp(mod(Ep,5),load_column)+ (fuel_mass_exp((5-mod(Ep,5)),load_column)-fuel_mass_exp(mod(Ep,5),load_column))*(Ep-(Ep-mod(Ep,5)))/((Ep+(5-mod(Ep,5)))-(Ep-mod(Ep,5)));          % fuel mass per cycle
end


mmix = (1 + AF_stoi)*mf;                        % mix (fuel + air) mass per cycle
%% Start and end CA per step

% 1. Intake
Ca_int_start    = 0;
Ca_int_end      = 180;

% 2. Compression
Ca_comp_start   = Ca_int_end;
Ca_comp_end     = 340;%340, also adjust in line 30 of dQwall_loss.m

%E0 NL 340 - 450, power = 0.53 kW
%E0 HL 340 - 425, power = 1.20 kW
%E0 FL 340 - 415, power = 2.15 kW

%E5 NL 340 - 435, power = 0.52 kW
%E5 HL 340 - 415, power = 1.17 kW
%E5 FL 340 - 400, power = 1.80 kW

%E10 NL 340 - 425
%E10 HL 340 - 410
%E10 FL 340 - 387

%E15 NL 340 - 425
%E15 HL 340 - 400
%E15 FL 340 - 385

% 3. Combustion
Ca_comb_start   = Ca_comp_end;      
Ca_comb_end     = 400;%400, also adjust in line 32 of dQwall_loss.m           
dCa_comb        = Ca_comb_end - Ca_comb_start;

% Note that the combustion ends before the exhaust valves open

% 4. Exhaust
Ca_ex_start     = 540;
Ca_ex_end       = 720;

%% Combustion heat release (Wiebe)

% Wiebe variables
a = 5;
n = 3;

den_g = 0.74;       % [kg/L] (only used in ratio, so unit not very important)
den_e = 0.78945;    % [kg/L]
den_mix = (den_g*(100-percentage) + den_e*percentage)./100;

mass_perc_E = percentage.*(den_e./den_mix);
Q_LHV = ((100 - mass_perc_E)./100).*43.4e3 + (mass_perc_E./100).*26.95e3;             % [kJ/kg]

theta_d = dCa_comb;         % ca difference start and end combustion 400-340
theta_s = Ca_comb_start;    % ca at start of combustion

% Defining angle range
theta = Ca_comb_start:1:Ca_comb_end;        % Wiebe function only applies in this region
theta(1) = theta_s;

NTheta = 360;               % Number of crank-angles
dTheta = 0.5;               % Stepsize
NSteps = NTheta/dTheta;

% Angle difference cannot be negative
if theta < theta_s;
   theta - theta_s  == 0;
end

for i = 1:length(theta)
    xb(i) = 1 - exp(-a*((theta(i)-theta_s)/theta_d)^n);                           % Fraction of burned fuel
    dQcom(i) = Q_LHV*mf*n*a*(1-xb(i))/theta_d*((theta(i)-theta_s)/theta_d)^(n-1); % Heat Release by combustion
    dQcom_T = transpose(dQcom);
end;

%% Real cycle
Ca=0:1:719;                          % defining CA length

V=zeros(length(Ca));                 % empty volume array with length CA
T=zeros(length(Ca));                 % empty temperature array with length CA
p=zeros(length(Ca));                 % empty pressure array with length CA

dQl=zeros(length(Ca));               % empty heat loss array with length CA
dQtot = 0;                           % adiabatic 

for dCa=1:length(Ca)
    V(dCa)= pi*(bore/2)^2*(rod+r_crank-(r_crank*cosd(dCa)+sqrt(rod^2-r_crank^2*(sind(dCa))^2)))+ v_c;          % Volume dependent on crank-angle
    
    % Defining a loop for each step to compute dQ, and subsequently T and p
    
    % Intake
    if (Ca_int_start <= dCa)&& (dCa <= Ca_int_end)   
        T(dCa) = T1;
        p(dCa) = p_int;
        
        dQwall(dCa) = dQwall_loss(dCa-1,dCa,T(dCa),p(dCa),p_int, gamma);
    end
    
    % Compression
 if (Ca_comp_start < dCa)&&(dCa < Ca_comp_end)                 
        [Cv_before, Cp_before, Rg_before] = before_comb(percentage/100, T(dCa-1));
        cv = Cv_before;
        R_g = Rg_before;
        gamma = Cp_before/Cv_before;

        dV = V(dCa)-V(dCa-1);                                   
        dQwall(dCa) = dQwall_loss(dCa-1,dCa,T(dCa-1),p(dCa-1),p_int, gamma);
        dU = dQwall(dCa) - ((p(dCa-1).*dV)./1000);                     % Via first law of thermodynamics
        dT = dU./(cv.*mmix);                                           
        
        T(dCa) = T(dCa-1) + dT;
        p(dCa) = (mmix.*R_g.*T(dCa))./V(dCa);
    end

    % Combustion
    if (Ca_comb_start <= dCa) && (dCa < Ca_comb_end)                
       [Cv_during, Cp_during, Rg_during] = during_comb(percentage/100, T(dCa-1), xb(dCa-Ca_comb_start+1));

        cv = Cv_during;
        R_g = Rg_during;
        gamma = Cp_during/Cv_during;
        
        dV = V(dCa)-V(dCa-1); 
        dQwall(dCa) = dQwall_loss(dCa-1,dCa,T(dCa-1),p(dCa-1),p_int, gamma);
        dU = dQcom(dCa-Ca_comb_start+1) +dQwall(dCa) - ...       % look at transpose wiebe
            ((p(dCa-1).*dV)./1000); 
     
        dT = dU./(cv.*mmix);                                           
        
        T(dCa) = T(dCa-1) + dT;
        p(dCa) = (mmix.*R_g.*T(dCa))./V(dCa);
    end  
     
    % Between end combustion and opening valves
    if (Ca_comb_end <= dCa) && (dCa < Ca_ex_start)       
       [Cv_after, Cp_after, Rg_after] = after_comb(percentage/100, T(dCa-1));

        cv = Cv_after;
        R_g = Rg_after;
        gamma = Cp_after/Cv_after;
        
        dV = V(dCa)-V(dCa-1);                                       
        dQwall(dCa) = dQwall_loss(dCa-1,dCa,T(dCa-1),p(dCa-1),p_int, gamma);
        dU = dQwall(dCa) - ((p(dCa-1).*dV)./1000);                    
        dT = dU./(cv.*mmix);                                           
        
        T(dCa) = T(dCa-1) + dT;
        p(dCa) = (mmix.*R_g.*T(dCa))./V(dCa);
    end
    
    % Exhaust valves open
    if (dCa == Ca_ex_start)   
        p(dCa) = Pamb;
        T(dCa) = (Pamb.*T(dCa-1))./p(dCa-1);
    end
    
    % Between closing exhaust valves and opening intake valves
    if (Ca_ex_start <= dCa) && (dCa < Ca_ex_end)
        p(dCa) = Pamb;
        T(dCa) = T(dCa-1);
    end
    
    % Intake valves open
    if (dCa == Ca_ex_end)                        
        p(dCa) = p_int;
        T(dCa)=(p_int.*T(dCa-1))./p(dCa-1);
    end
    
end

figure()
plot(V,p*10^-5,'LineWidth',6)
set(gca,'FontSize',30)
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
title({sprintf('PV Diagram, E%.0f %s load', [percentage, load])}, 'FontSize', 50)
grid on
give_pv_diagram(percentage, load);
if ismember(percentage,[0, 5, 10, 15])
    legend("Model","Experiment")
end
%%
angle = 0:1:719;
figure()
plot(angle,p*10^-5,'LineWidth',6)
set(gca,'FontSize',30)
xlabel('Crank angle [\theta]')
ylabel('Pressure [Bar]')
xlim([0 720])
title({sprintf('Pressure vs crank angle, E%.0f %s load', [percentage, load])}, 'FontSize', 50)
grid on

%%
%Work
p_V_array = [V(:,1) p(:,1)];
p_V_array(p_V_array(:, 2) > Pamb, :) = [];

total_work = polyarea(V,p);
total_work = total_work(1);
negative_work = polyarea(p_V_array(:,1),p_V_array(:,2));

work = total_work - 2*negative_work; %[J]


measured_voltage = 203;
%Determines the amperes. No load is 0, half load is 4 and full load is 7.3 (See Arjan&Sven SSA)
if load == "no"
    P_out = measured_voltage * 0;
elseif load == "half"
    P_out = measured_voltage * 4;
elseif load == "full"
    P_out = measured_voltage * 7.3;
end
V_dis = pi/4 * bore^2 * stroke * 4;
p_me = (p1+p2+p3+p4)/4;
W = p_me * V_dis;
P_in = W * (50/2);
mechanical_efficiency = P_out / P_in;


% cv = cv_function(percentage,T)%cv based on percentage ethanol and temperature
% T4=(p4/p1)*T1%Calculated with ideal gas law
% Qc = cv*(MAir+MFuel)*(T4-T1)
% 
% %efficiency
% efficiency = work/(work+Qc)
% return
%%
%Power
dt_no = 0.0378;%Time per cycle, no load
dt_half = 0.0389344545;%Time per cycle, half load
dt_full = 0.040065;%Time per cycle, full load

if load == "no"
    power = work/dt_no; %J/s is equal to W. Power per cycle
elseif load == "half"
    power = work/dt_half; %J/s is equal to W. Power per cycle
elseif load == "full"
    power = work/dt_full; %J/s is equal to W. Power per cycle
end

disp((sprintf('%.2f kW', [power/1000])))  % kW