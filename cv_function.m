function [cv_mix] = cv_function(percentage,T)

%% To make sure that matlab will find the functions. You must change it to your situation 
abspath_to_generalfolder='General'; % absolute reference to General folder
addpath(abspath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase); %#ok<LOAD>
global Runiv
Runiv = 8.314;

%% 
DBdir = 'General\Nasa';
DBname = 'NasaThermalDatabase';
load(fullfile(DBdir,DBname));
%% Add ethanol to database
iCHO=myfind({El.Name},{'C','H','O'}); % Which position corresponds to C H and O elements in the Sp.Elcomp vector
iH=iCHO(2);iC=iCHO(1);iO=iCHO(3);
iAlc=[];ii=0;
for i=1:length(Sp)
    if ( Sp(i).Elcomp(iH)==(2*Sp(i).Elcomp(iC)+2) && Sp(i).Elcomp(iO)==1) % Typical comp of an alcohol 2*nC+2 H-atoms, 1 O-atom
        ii = ii+1;
        iAlc(ii)=i;
    end;
end
%% Find index of species
iSp = myfind({Sp.Name},{'Gasoline','O2','N2','CO2','H2O','C2H5OH'});

% Cv
cv = [CvNasa(T,Sp(iSp(1))), CvNasa(T,Sp(iSp(2))), CvNasa(T,Sp(iSp(3))), CvNasa(T,Sp(iSp(4))), CvNasa(T,Sp(iSp(5))), CvNasa(T,Sp(iSp(6)))];

%% Stoichiometric combustion
% (1-percentage)*Cx_gasHy_gasOz_gas + percentage*Cx_ethHy_ethOz_eth +
% ((x_gas+x_eth)+(y_gas+y_eth)/4-(z_gas+z_eth)/2)(O2 + 3.76N2) -->
% (x_gas+x_eth)CO2 + ((y_gas+y_eth)/2)H2O + ((x_gas+x_eth)+(y_gas+y_eth)/4-(z_gas+z_eth)/2)N2

% mol_gas = (100-percentage)*volume*density * molaire massa
mol_gas = (100-percentage)*1*0.75*/(Sp(iSp(1)).Mass*1000);                        % Density gasoline is from 0.71 to 0.77 g/cm3

% mol_ethanol = percentage*volume*density * molaire massa
mol_ethanol = percentage*1*0.78945/(46);                                     % Density Ethanol    0.78945 g/cm3 (at 20 ?C)

x_gas = mol_gas*Sp(iSp(1)).Elcomp(3); 
y_gas = mol_gas*Sp(iSp(1)).Elcomp(2); 
z_gas = mol_gas*Sp(iSp(1)).Elcomp(1); 

x_eth = mol_ethanol*2;
y_eth = mol_ethanol*6;
z_eth = mol_ethanol*1;

%% Molair masses
C = Sp(myfind({Sp.Name},{'C'})).Mass * 1000;      % Carbon    [g/mol]
H = Sp(myfind({Sp.Name},{'H'})).Mass * 1000;      % Hydrogen  [g/mol]
O = Sp(myfind({Sp.Name},{'O'})).Mass * 1000;      % Oxygen    [g/mol]
N = Sp(myfind({Sp.Name},{'N'})).Mass * 1000;      % Nitrogen  [g/mol]
%% Masses
M_gas = C*x_gas + H*y_gas + O*z_gas;                        % [g]
M_eth = C*x_eth + H*y_eth + O*z_eth;                        % [g]

M_gas_air = (x_gas+y_gas/4-z_gas/2)*(O*2 + 3.76*N*2);   % [g]
MO_gas = (x_gas+y_gas/4-z_gas/2)*O*2;                  % [g]

M_eth_air = (x_eth+y_eth/4-z_eth/2)*(O*2 + 3.76*N*2);   % [g]
MO_eth = (x_eth+y_eth/4-z_eth/2)*O*2;                  % [g]

MN_gas = 3.76*(x_gas+(y_gas/4)-(z_gas/2))*N*2;             % [g]
MN_eth = 3.76*(x_eth+(y_eth/4)-(z_eth/2))*N*2;             % [g]

%% Cv for mixtures
M_mix_gas = M_gas + M_gas_air;  % total mass reactants 0% ethanol (100% gasoline) combustion
M_mix_eth = M_eth + M_eth_air;  % total mass reactants 100% ethanol combustion

% cv for gasoline and ethanol mixtures [J/KgK]
cv_mix_gas = (M_gas/M_mix_gas)*cv(1)+(MO_gas/M_mix_gas)*cv(2)+(MN_gas/M_mix_gas)*cv(3);
cv_mix_eth = min(0,(M_eth/M_mix_eth)*cv(6)+(MO_eth/M_mix_eth)*cv(2)+(MN_eth/M_mix_eth)*cv(3));

% cv for blends of ethanol and gasoline [kJ/kgK]
cv_mix = ((100-percentage)*cv_mix_gas + percentage*cv_mix_eth)./1000;

