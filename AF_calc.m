%% Stoichiometric combustion
% CxHyOz + (x + y/4-z/2)(O2) --> xCO2 + (y/2)H2O
x_gas = 7.76;
y_gas = 13.1;
z_gas = 0;

x_eth = 2;
y_eth = 6;
z_eth = 1;
%% Molair masses
C = 12.01;      % Carbon    [g/mol]
H = 1.008;      % Hydrogen  [g/mol]
O = 16.00;      % Oxygen    [g/mol]
N = 14.01;      % Nitrogen  [g/mol]
%% AF Gasoline
M_gas = C*x_gas + H*y_gas + O*z_gas;                        % [g]
M_air_gas = (x_gas+(y_gas/4)-(z_gas/2))*(O*2 + 3.76*N*2);   % [g]
MO2_gas = (x_gas+(y_gas/4)-(z_gas/2))*O*2;                  % [g]
AF_st_gas = M_air_gas/M_gas

M_eth = C*x_eth + H*y_eth + O*z_eth;                        % [g]
M_air_eth = (x_eth+(y_eth/4)-(z_eth/2))*(O*2 + 3.76*N*2);   % [g]
MO2_eth = (x_eth+(y_eth/4)-(z_eth/2))*O*2;                  % [g]
AF_st_eth = M_air_eth/M_eth