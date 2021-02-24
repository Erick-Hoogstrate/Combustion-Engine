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
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel

density_gasoline = 1.319;                                                    %kg/liter
mfurate = 1.3*density_gasoline*1/3600*kg/s;  
r = 8.5;

x= SpS(1).Elcomp(3);                                                        %The amount of moles of carbon in the reactants, equal to the amount of moles of C02 in the products
y= SpS(1).Elcomp(2);                                                        %the amount of moles of hydrogen in the reactants                                             
O2_consumed = (x + (y/4))                                                   %The amount of moles of 02 in the air consumed by the ideal stochiometric combustion proccess,from eq 3.29 in Turns
MFuel = SpS(1).Mass;                                                        %defining the molar mass of the fuel in a variable from the SpS struct
AF_stoic = 4.76 * O2_consumed * ( MAir / MFuel);                            %dirrectly from equation 3.30 in the book
AF = AF_stoic;
%%
%state variables
p1 = Pamb;
T1 = Tamb;

%input
T3=2300;                                                                  %unclear
%T3 = 503;
gamma = 1.4;

%engine geometric parameters
bore = 0.068;
stroke = 0.054;
rod = 0.091313;
r = 8.5;
RPM = 3000;

%calculating swept volume and the clearance volume
v_d = (pi/4)*bore^2*stroke;
v_c = v_d/(r-1); 

%starting point
p0 = p1;
v0 = (v_d + v_c)/r

%point 1
v1 = v_d + v_c;
v2 = v_c;
g = 1.4;    %gamma

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
P0=plot([v0 v1],[p0 p1],'-','color','blue');%isobaric expansion
P1=plot(v_compression,p_compression,'color','m');%isentropic compression
P2=plot([v2 v3],[p2 p3],'color','b');%constant volume heat addition
P3=plot(v_expansion,p_expansion,'color','k');%isentropic expansion
P4=plot([v4 v1],[p4,p1],'color','g');%constant volume heat rejection
P5=plot([v1 v0],[p1 p0],'--','color','yellow');%isobaric compresion
xlabel('volume ')
ylabel('pressure')
title({'PV Diagram';'Otto cycle'})
legend([P0 P1 P2 P3 P4 P5],{'isobaric expansion','isentropic compression','constant volume heat addition','isentropic expansion','constant volume heat rejection', 'isobaric compression'})

%thermal efficiency
thermalefficiency=(1-1/r^(g-1))    

%%
%Work
area12 = trapz(-v_compression,p_compression)                                 %start for calculation work
area34 = trapz(v_expansion,p_expansion)

work = area34-area12%Work per cycle

%%
%Heat losses

%constants
Tw = 323;      % Estimation of the wall temperature of the engine (50 degrees celcius)
C0 = 120;    % between 110-130
% P =       % dependent on crank angle function
% A =       % dependent on crank angle function
V_mp = 2*stroke*(RPM/60);

% h_g = C0*bore^-0.2*P^0.8*((C1*V_mp)+(C2*(V_d*T1)/(p1*v1)*(P-P_mot)))^0.8*T^-0.53

% step 0-1 (isobaric expansion)
C1 = 6.18;
C2 = 0;
P = p0;
T = T1;

h_g_01 = C0*bore^-0.2*P^0.8*((C1*V_mp))^0.8*T^-0.53

% dQdt = -h_g_01*A*(T-Tw)


% step 1-2 (isentropic compression)
C1 = 2.28;
C2 = 0;
P = p_compression;
T = T2;
% T = (p_compression.*v_compression*T1/(p1*v1);     % Unclear what T should be used and whether this should be a constant value or a function

h_g_12 = C0*bore^-0.2*P.^0.8.*((C1*V_mp))^0.8*T^-0.53


% dQdt = -h_g*A*(T3-Tw)



return
cv=cp/gamma_cv%Use NasaUseExample.m to retrieve cp and gamma. From this you can get cv. Gamma specific to substance.
T4=(p4/p1)*T1%Calculated with ideal gas law
Qc = cv*(MAir+MFuel)*(T4-T1)



%%
%efficiency

efficiency = work/(work+Qc)






return
%%
%Power
dt_no = 0.0378%Time per cycle, no load
dt_full = 0.040065%Time per cycle, full load
dt_half = 0.0389344545%Time per cycle, half load

power = work/dt_no%Power per cycle

return
%%
% Initialisation
p(1)=P0;T(1)=T0;
pad(1)=p(1);
Ca(1)=0.0;V(1)=Vcyl(Ca(1),S,B,l,rc); % Vcyl is a function that
% computes cyl vol as fie of crank-
% angle for given B,S,l and rc
m(1)=p(1)*V(1)/Rg/T(1);
%
%
% Loop over crank-angle, with 'for' construction
NCa=360; % Number of crank-angles
dCa=0.5; % Stepsize
NSteps=NCa/dCa;
for i=2:NSteps,
    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc); % New volume for current crank-angle
    m(i)=m(i-1); % Mass is constant, valves are closed
    dV=V(i)-V(i-1); % Volume change
    dQcom = YourModel(Ca(i)); % Heat Release by combustion
    dT=(dQcom-p(i-1)*dV)/Cv/m(i-1); % 1st Law dU=dQ-pdV (closed system)
    % adiabatic closed system with constant
    % gas composition and constant Cv
    T(i)=T(i-1)+dT;
    p(i)=m(i)*Rg*T(i)/V(i); % Gaslaw
end;
