function [Cv_before, Cp_before, Rg_before] = before_comb(E, Tr)
abspath_to_generalfolder='Nasa'; % absolute reference to General folder
addpath(abspath_to_generalfolder); 

TdataBase=fullfile('General', 'Nasa','NasaThermalDatabase');
load(TdataBase); 
global Runiv
Runiv = 8.314;
%% 
iSp = myfind({Sp.Name},{'Gasoline','O2','N2','CO2','H2O','C2H5OH'});

C = 12.01;     
H = 1.008;     
O = 16.00;      
N = 14.01;      

CO2 = C+2*O;    
H2O = 2*H + O;  
O2  = 2*O;       
N2  = 2*N; 

T = Tr(1);

Cv = [CvNasa(T,Sp(iSp(1))), CvNasa(T,Sp(iSp(2))), CvNasa(T,Sp(iSp(3))), CvNasa(T,Sp(iSp(4))), CvNasa(T,Sp(iSp(5))), CvNasa(T,Sp(iSp(6)))];
Cp = [CpNasa(T,Sp(iSp(1))), CpNasa(T,Sp(iSp(2))), CpNasa(T,Sp(iSp(3))), CpNasa(T,Sp(iSp(4))), CpNasa(T,Sp(iSp(5))), CpNasa(T,Sp(iSp(6)))];
%% 
x_gas = 7.76;
y_gas = 13.1;
z_gas = 0;

x_eth = 2;
y_eth = 6;
z_eth = 1;

M_gas     = C*x_gas + H*y_gas + O*z_gas;                        
M_air_gas = (x_gas+(y_gas/4)-(z_gas/2))*(O2 + 3.76*N2);   
M_O2_gas  = (x_gas+(y_gas/4)-(z_gas/2))*O2;                
AF_st_gas = M_air_gas/M_gas;

M_eth     = C*x_eth + H*y_eth + O*z_eth;                     
M_air_eth = (x_eth+(y_eth/4)-(z_eth/2))*(O2 + 3.76*N2);   
M_O2_eth   = (x_eth+(y_eth/4)-(z_eth/2))*O2;                  
AF_st_eth = M_air_eth/M_eth;

M_CO2_gas = x_gas*(CO2);                         % [g]
M_H2O_gas = (y_gas/2)*(H2O);                     % [g]
M_N2_gas = 3.76*(x_gas+(y_gas/4)-(z_gas/2))*N2;  % [g]

M_CO2_eth = x_eth*(CO2);                         % [g]
M_H2O_eth = (y_eth/2)*(H2O);                     % [g]
M_N2_eth = 3.76*(x_eth+(y_eth/4)-(z_eth/2))*N2;  % [g]
%%
M_mix_gas = M_gas + M_air_gas;  
M_mix_eth = M_eth + M_air_eth;  

Cv_mix_gas = (M_gas/M_mix_gas)*Cv(1)+(M_O2_gas/M_mix_gas)*Cv(2)+(M_N2_gas/M_mix_gas)*Cv(3);
Cp_mix_gas = (M_gas/M_mix_gas)*Cp(1)+(M_O2_gas/M_mix_gas)*Cp(2)+(M_N2_gas/M_mix_gas)*Cp(3);
Cv_mix_eth = (M_eth/M_mix_eth)*Cv(6)+(M_O2_eth/M_mix_eth)*Cv(2)+(M_N2_eth/M_mix_eth)*Cv(3);
Cp_mix_eth = (M_eth/M_mix_eth)*Cp(6)+(M_O2_eth/M_mix_eth)*Cp(2)+(M_N2_eth/M_mix_eth)*Cp(3);

Cv_before = ((1-E)*Cv_mix_gas + E*Cv_mix_eth)./1000;   
Cp_before = ((1-E)*Cp_mix_gas + E*Cp_mix_eth)./1000;   

Rg_before = (Cp_before-Cv_before)*1000;       