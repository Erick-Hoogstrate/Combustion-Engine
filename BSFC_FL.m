%BSFC = fuel consumption (g/s) / Power produced (W)
%Output: 7.33 Ampère & 212 volt
%% E0 BSFC FL
dens_E0 = 0.74; %[g/cm^3]
%takes 8.20s to burn 4ml so 8.20s to burn 4cm3
time = 8.2;
BSFCE0=(0.74/(8.2/4))/(212*7.33);

%% E5 BSFC FL
dens_E5 = (0.74*95 + 0.78945*5)/100;
mass_E5 = dens_E5* 10^-3 * 4;
time = 9.67;
BSFCE5=(dens_E5/(9.67/4))/(212*7.33);

%% E10 BSFC FL
dens_E10 = (0.74*90 + 0.78945*10)/100;
mass_E10 = dens_E10* 10^-3 * 4;
time = 10.82;
BSFCE10=(dens_E10/(10.82/4))/(212*7.33);

%% E15 BSFC FL
dens_E15 = (0.74*85 + 0.78945*15)/100;
mass_E15 = dens_E15* 10^-3 * 4;
time = 11.14;
BSFCE15=(dens_E15/(11.14/4))/(212*7.33);

%% Turn into [g/kWh]
%Answers were now in [g/J]
%turn to g/kWh

BSFCE0*(3.6*10^6)   %836.2585
BSFCE5*(3.6*10^6)   %711.5027
BSFCE10*(3.6*10^6)  %637.9984
BSFCE15*(3.6*10^6)  %621.7284

%%
%eff=1/(BSFC*LHV in kwh/g)
%LHV gasoline=43.44 MJ/kg=12.0667 kwh/kg=0.0120667 kwh/g
% LHVgas = 0.0120667;
% LHVeth = 0.0074861;
% eff5=1/(BSFCE5*LHVeth) % = 6.7588e5




