%% E0
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [21201, 25200];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["E0", "E1", "E2"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
TestE0noload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE0_no_load_2.txt", opts);

t_nl_0            = TestE0noload2.E0;
puls_sens_nl_0  = TestE0noload2.E1;
pres_sens_nl_0 = TestE0noload2.E2;

%% E5
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [21001, 25000];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["E0", "E1", "E2"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

TestE5noload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE5_no_load_2.txt", opts);

t_nl_5            = TestE5noload2.E0;
puls_sens_nl_5  = TestE5noload2.E1;
pres_sens_nl_5 = TestE5noload2.E2;

%% E10
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [18501, 22500];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["E0", "E1", "E2"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
TestE10noload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE10_no_load_2.txt", opts);

t_nl_10            = TestE10noload2.E0;
puls_sens_nl_10  = TestE10noload2.E1;
pres_sens_nl_10 = TestE10noload2.E2;

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
TestE15noload2 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\University\4GB10\raw_data_exp2\TestE15_no_load_2.txt", opts);

t_nl_15            = TestE15noload2.E0;
puls_sens_nl_15  = TestE15noload2.E1;
pres_sens_nl_15 = TestE15noload2.E2;
%%
pressure_relative_nl_0  = (pres_sens_nl_0 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_nl_0 = pressure_relative_nl_0 - 1.165  % ?? estimated (relative sensor)

pressure_relative_nl_5  = (pres_sens_nl_5 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_nl_5 = pressure_relative_nl_5 % ?? estimated (relative sensor)

pressure_relative_nl_10  = (pres_sens_nl_10 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_nl_10 = pressure_relative_nl_10 - 1.6   % ?? estimated (relative sensor)

pressure_relative_nl_15  = (pres_sens_nl_15 -(0.115*5))/(0.00385*5*4); %formula from Canvas 
pressure_nl_15 = pressure_relative_nl_15  % ?? estimated (relative sensor)

phi_speed = 18000*t_nl_0
phi = phi_speed 

%% One figure has all
figure()
plot(phi, pressure_nl_0)
hold on
plot(phi, pressure_nl_5)
hold on
plot(phi, pressure_nl_10)
hold on
plot(phi, pressure_nl_15)
legend('E0', 'E5', 'E10', 'E15')
xlabel('Crank Angle [phi]')
ylabel('Pressure [bar]')
title('No Load')

%%
a = 0.027; %radius (equal to stroke/2)
B = 0.068; %diameter of bore
L = 0.085; %length of rod
Vc = 5/1000000; %deathvolume or clearance volume, estimated value %aanpassen

s = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); %explained in ssa
Volume = (((pi)*B.^2)/(4))*(L + a - s)+ Vc; %explained in ssa

%% One figure has all pv
figure()
plot(Volume, pressure_nl_0)
hold on
plot(Volume, pressure_nl_5)
hold on
plot(Volume, pressure_nl_10)
hold on
plot(Volume, pressure_nl_15)
legend('E0', 'E5', 'E10', 'E15')
xlabel('Volume [m^3]')
ylabel('Pressure [bar]')
title('No Load PV')

%% Multiple figures
% Create plots
% 
% figure('Name','No Load')
% subplot(1,4,1)
% plot(phi, pressure_nl_0)
% title('E0 no load')
% ylim([0,10])
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')
% subplot(1,4,2)
% plot(phi, pressure_nl_5)
% title('E5 no load')
% ylim([0,10])
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')
% subplot(1,4,3)
% plot(phi, pressure_nl_10)
% title('E10 no load')
% ylim([0,10])
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')
% subplot(1,4,4)
% plot(phi, pressure_nl_15)
% title('E15 no load')
% ylim([0,10])
% xlabel('Crank Angle')
% ylabel('Pressure (bar)')

%PV 

% figure('Name','No Load')
% subplot(2,2,1)
% plot(Volume, pressure_nl_0)
% title('E0 no load')
% ylim([0,8])
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure (bar)')
% subplot(2,2,2)
% plot(Volume, pressure_nl_5)
% title('E5 no load')
% ylim([0,10])
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure (bar)')
% subplot(2,2,3)
% plot(Volume, pressure_nl_10)
% title('E10 no load')
% ylim([0,10])
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure (bar)')
% subplot(2,2,4)
% plot(Volume, pressure_nl_15)
% title('E15 no load')
% ylim([0,10])
% xlim([0,2.1e-4])
% xlabel('Volume')
% ylabel('Pressure (bar)')

%% Clear temporary variables
clear opts