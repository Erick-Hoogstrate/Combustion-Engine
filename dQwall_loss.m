function [dQ_loss] = dQwall_loss(Ca1,Ca2,T,p,p_int,gamma)

% Constants
Tamb = 298.15;              % Temperature enviroment
Tw   = 363;                 % Wall temperature, assumed 90 degrees
stroke = 0.054;             % [m]
rod  = 0.091313;            % [m]
bore = 0.068;               % [m]

r_crank = stroke/2;         % Radius of the crankshaft[m] (half the stroke)
r_cyl   = bore/2;           % Radius cyllinder [m]

RPM = 3000;
MPS = 2*stroke*RPM/60;      % Mean piston speed [m/s]
gamma = 1.4;                % Gamma
r = 8.5;                    % Compression ratio


% Volume
v_d = (pi/4)*bore^2*stroke; % volume BDC
v_c = v_d/(r-1);            % clearance volume

V = pi*(bore/2)^2*(rod+r_crank-(r_crank*cosd(Ca1)+sqrt(rod^2-r_crank^2*(sind(Ca1))^2))) + v_c;      % Volume dependent on CA

% Woschni heat transfer coefficients depend on CA

% These values are the same as main model
Ca_int_end      = 180;
Ca_comp_start   = Ca_int_end;
Ca_comp_end     = 340;
Ca_comb_start   = Ca_comp_end;      
Ca_comb_end     = 400; 

% Coefficients and formulas have to be validated

if (Ca_comb_start<=Ca2) && (Ca2<Ca_comb_end)
    C1=2.28;
    C2=3.24*10^-3;
   
elseif (Ca_comp_start<Ca2)&&(Ca2<Ca_comp_end)
    C1=2.28;
    C2=0;
    
elseif (Ca2<=Ca_int_end) || (Ca2>=Ca_comb_end)
    C1=6.18;
    C2=0;
end

V_ref = v_d;                % reference volume equals V at BDC
P_amb=1.01235e5;
p_ref = P_amb;              % ambient pressure

% Heat transfer coefficient calculation
p_m = p_ref.*(V_ref./V).^gamma;
v = C1.*MPS+C2.*((v_d.*Tamb)./(p_ref.*V_ref)).*(p-p_m);
coefficient = 0.013.*(bore^-0.2).*(p.^0.8).*(T.^-0.53).*(v.^0.8);

A_1 = pi*r_cyl^2; 
A_2 = 2*pi*r_cyl*(rod+r_crank-(r_crank*cosd(Ca1)+sqrt(rod^2-r_crank^2*(sind(Ca1))^2))); 
A_tot = 2*A_1 + A_2; 

% Heat loss 
dQ = -coefficient.*(T-Tw).*A_tot; 
dt = (1/(RPM/60))/360;                      % duration per CA
dQ_loss = dQ.*dt.*(Ca2-Ca1)./1000;          % total heat loss [kJ]