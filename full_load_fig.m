%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\20192420\Documents\Q3 Y2\data exp 2\TestE10_full_load_2.txt
%
% Auto-generated by MATLAB on 03-Mar-2021 09:49:14

%% Setup the Import Options (E0)
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [18451, 22450];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["E0", "E1", "E2"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
TestE0fullload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE0_full_load_2.txt", opts);

t_fl_0            = TestE0fullload2.E0;
puls_sens_fl_0  = TestE0fullload2.E1;
pres_sens_fl_0 = TestE0fullload2.E2;

%% E5
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [20251, 24250];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["E0", "E1", "E2"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

TestE5noload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE5_full_load_2.txt", opts);

t_fl_5            = TestE5noload2.E0;
puls_sens_fl_5  = TestE5noload2.E1;
pres_sens_fl_5 = TestE5noload2.E2;

%% E10
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [20101, 24100];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["E0", "E1", "E2"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
TestE10noload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE10_full_load_2.txt", opts);

t_fl_10            = TestE10noload2.E0;
puls_sens_fl_10  = TestE10noload2.E1;
pres_sens_fl_10 = TestE10noload2.E2;

%% E15
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [19001, 23000];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["E0", "E1", "E2"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
TestE15noload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE15_full_load_2.txt", opts);

t_fl_15            = TestE15noload2.E0;
puls_sens_fl_15  = TestE15noload2.E1;
pres_sens_fl_15 = TestE15noload2.E2;
%%
pressure_relative_fl_0  = (pres_sens_fl_0 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_fl_0 = pressure_relative_fl_0 -23   % ?? estimated (relative sensor)

pressure_relative_fl_5  = (pres_sens_fl_5 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_fl_5 = pressure_relative_fl_5 -8.4 % ?? estimated (relative sensor)

pressure_relative_fl_10  = (pres_sens_fl_10 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_fl_10 = pressure_relative_fl_10 -18.7    % ?? estimated (relative sensor)

pressure_relative_fl_15  = (pres_sens_fl_15 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_fl_15 = pressure_relative_fl_15 - 15.8 % ?? estimated (relative sensor)

phi_speed = 18000*t_fl_0 
phi = phi_speed %

% figure()
% plot(phi, pressure_fl)

%% Multiple figures Create plots

% figure('Name','Full load')
% subplot(1,4,1)
% plot(phi, pressure_fl_0)
% title('E0 full load')
% 
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')
% subplot(1,4,2)
% plot(phi, pressure_fl_5)
% title('E5 full load')
% 
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')
% subplot(1,4,3)
% plot(phi, pressure_fl_10)
% title('E10 full load')
% 
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')
% subplot(1,4,4)
% plot(phi, pressure_fl_15)
% title('E15 full load')
% 
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')

%% One figure has all

figure()
plot(phi, pressure_fl_0)
hold on
plot(phi, pressure_fl_5)
hold on
plot(phi, pressure_fl_10)
hold on
plot(phi, pressure_fl_15)
legend('E0', 'E5', 'E10', 'E15')
xlabel('Crank Angle [phi]')
ylabel('Pressure [bar]')
title('Full Load')
%%
%

a = 0.027; %radius (equal to stroke/2)
B = 0.068; %diameter of bore
L = 0.085; %length of rod
Vc = 5/1000000; %deathvolume or clearance volume, estimated value %aanpassen

s = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); %explained in ssa
Volume = (((pi)*B.^2)/(4))*(L + a - s)+ Vc; %explained in ssa

% figure('Name','full load')
% subplot(2,2,1)
% plot(Volume, pressure_fl_0)
% title('E0 full load')
% 
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure (bar)')
% subplot(2,2,2)
% plot(Volume, pressure_fl_5)
% title('E5 full load')
% 
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure (bar)')
% subplot(2,2,3)
% plot(Volume, pressure_fl_10)
% title('E10 full load')
% 
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure [bar]')
% subplot(2,2,4)
% plot(Volume, pressure_fl_15)
% title('E15 full load')
% 
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure [bar]')

%% One figure has all pv

figure()
plot(Volume, pressure_fl_0)
hold on
plot(Volume, pressure_fl_5)
hold on
plot(Volume, pressure_fl_10)
hold on
plot(Volume, pressure_fl_15)
legend('E0', 'E5', 'E10', 'E15')
xlabel('Volume [m^3]')
ylabel('Pressure [bar]')
title('Full Load PV')
% 
% % Clear temporary variables
clear opts