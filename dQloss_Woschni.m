function [dQ] = dQlossW(T,theta1,theta2,p,p_intake,g)
% This function is only to be used with no more than 1 degrees crankangle difference!
Tamb=298.15; % Reference Temperature
Tw = 273+90; % wall temperature = 90 degrees celcius
rod = 0.091313; %[m]
bore = 0.068; %[m]
radius = bore/2; %[m]
stroke = 0.054; %[m]
R = stroke/2; % Radius of the crankshaft[m] (half the stroke)
RPM = 3000;
MPS = 2*stroke*RPM/60;  % Mean piston speed [m/s]
r = 8.5; %compression ratio
g = 1.4; %gamma

%calculating swept volume and the clearance volume
%v_d =(pi/4)*bore^2*(rod+radius-(radius*cosd(180)+sqrt(rod^2-radius^2*(sind(180))^2))); %[m^3]
v_d = (pi/4)*bore^2*stroke; %[m^3]
v_c = v_d/(r-1); %[m^3]
%v_c = 25576.42*10^(-9); %[m^3]

V=v_c+(pi/4)*bore^2*(rod+R-(R*cosd(theta1)+sqrt(rod^2-R^2*(sind(theta1))^2)));  % Current volume combustion chamber

% Angles input ---- WOSCHNI
IntakeStop = 180;
CompressionStart = IntakeStop;
CompressionStop = 360;
CombustionStart = CompressionStop; 
CombustionStop = 430; %% Actual = 480, but literature says >60 degrees reaches >90% burnt

if (theta2<=IntakeStop) || (theta2>=CombustionStop)
    C1=6.18;
    C2=0;
elseif (CompressionStart<theta2)&&(theta2<CompressionStop)
    C1=2.28;
    C2=0;
elseif (CombustionStart<=theta2) && (theta2<CombustionStop)
    C1=2.28;
    C2=3.24*10^-3;
end

pref = p_intake;
Vref = v_c+(pi/4)*bore^2*(rod+R-(R*cosd(IntakeStop)+sqrt(rod^2-R^2*(sind(IntakeStop))^2)));  % Current volume combustion chamber

pm = pref.*(Vref./V).^g;
v = C1.*MPS+C2.*((v_d.*Tamb)./(pref.*Vref)).*(p-pm);
alpha = 0.01297829376.*(bore^-0.2).*(p.^0.8).*(T.^-0.53).*(v.^0.8)

A1 = pi*radius^2; % area of pistonhead / cylinderhead (inside), approximation
A2 = 2*pi*radius*(rod+R-(R*cosd(theta1)+sqrt(rod^2-R^2*(sind(theta1))^2))); % current area of cylinder wall
Atot = 2*A1+A2; % area pistondeck + area cylinderhead + area cylinder wall (approximation)

dQdt = -alpha.*(T-Tw).*Atot; % heat loss per second
dt = (1/(RPM/60))/360; % time taken per crankangle [s]
dQ = dQdt.*dt.*(theta2-theta1)./1000; % heat loss [kJ]