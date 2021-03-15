%BSFC = fuel consumption (g/s) / Power produced (W)
%Output: 4.06 Ampère & 224 volt
%% E0 BSFC HL
dens_E0 = 0.74; %[g/cm^3]
mass_E0 = dens_E0 *10^-3 * 4; %[g/cm^3]
%takes 8.20s to burn 4ml so 8.20s to burn 4cm3
time = 11.73;
BSFCE0=(0.74/(11.73/4))/(224*4.06)

%% E5 BSFC HL
dens_E5 = (0.74*95 + 0.78945*5)/100;
mass_E5 = dens_E5* 10^-3 * 4;
time = 12.59;
BSFCE5=(dens_E5/(12.59/4))/(224*4.06)

%% E10 BSFC HL
dens_E10 = (0.74*90 + 0.78945*10)/100;
mass_E10 = dens_E10* 10^-3 * 4;
time = 12.95;
BSFCE10=(dens_E10/(12.95/4))/(224*4.06)

%% E15 BSFC HL
dens_E15 = (0.74*85 + 0.78945*15)/100;
mass_E15 = dens_E15* 10^-3 * 4;
time = 14.74;
BSFCE15=(dens_E15/(14.74/4))/(224*4.06)

%% Turn into [g/kWh]
%Answers were now in [g/J]
%turn to g/kWh

BSFCE0*(3.6*10^6)   %998.9003
BSFCE5*(3.6*10^6)   %933.7768
BSFCE10*(3.6*10^6)  %910.8416
BSFCE15*(3.6*10^6)  %802.8866