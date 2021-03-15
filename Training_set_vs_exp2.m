%% -------------------------NO LOAD E0---------------------------------
clear all;close all;clc;

for i =0:9;
  
    opts = delimitedTextImportOptions("NumVariables", 3);
    
    rotation_speed = 48.75769465; %1 rotation = 0.0205095833s  
    samples        = fix((1/rotation_speed)*100000*2); 
    a              = samples * i;
    opts.DataLines = [1+a, samples+a]; 

    opts.Delimiter        = "\t";
    opts.VariableNames    = ["E0", "E1", "E2"];
    opts.VariableTypes    = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule    = "read";

    Test = readtable("D:\Documents\MATLAB\(4GB10) DBL Combustion engine\Experiment 2\raw_data_exp2\TestE0_no_load_3.txt", opts);

    t          = Test.E0;
    puls_sens  = Test.E1;
    pres_sens  = Test.E2;
    
    x = -1.4; %correcting pressure sensor
    
    while i == 0
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c1 = pressure_relative + x;
        i = i+1
    end
    
    while i == 1
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c2 = pressure_relative + x;
        i = i+1
    end
    
    while i == 2
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c3 = pressure_relative + x;
        i = i+1
    end
    
    while i == 3
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c4 = pressure_relative + x;
        i = i+1
    end
    
    while i == 4
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c5 = pressure_relative + x;
        i = i+1
    end
    
    while i == 5
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c6 = pressure_relative + x;
        i = i+1
    end
    
    while i == 6
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c7 = pressure_relative + x;
        i = i+1
    end
    
    while i == 7
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c8 = pressure_relative + x;
        i = i+1
    end
    
    while i == 8
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c9 = pressure_relative + x;
        i = i+1
    end
    
    while i == 9
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c10 = pressure_relative + x;
        i = i+1
    end
    
    clear opts
end

average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

start_angle  = 72.1419;
double_tooth = 83.13;
phi_speed    = rotation_speed*360*t; 
phi          = phi_speed-(start_angle+double_tooth);

%Training set no load
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines        = [1850, 5950];
opts.Delimiter        = "\t";
opts.VariableNames    = ["E0", "E1", "E2"];
opts.VariableTypes    = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule    = "read";

Training_set = readtable("D:\Documents\MATLAB\(4GB10) DBL Combustion engine\Test data\Training set\No_load.txt", opts);

t_relative    = Training_set.E0;
puls_sens_ts  = Training_set.E2;
pres_sens_ts  = Training_set.E1;

pressure_relative_training_set = (pres_sens_ts -(0.115*5))/(0.00385*5*4); %formula from Canvas
pressure_training_set = pressure_relative_training_set*4.5-0.481;
t_training_set = t_relative + 0.3508;

figure()
hold on 
% plot(t, p_c1)
% plot(t, p_c2)
% plot(t, p_c3)
% plot(t, p_c4) 
% plot(t, p_c5) 
% plot(t, p_c6)
% plot(t, p_c7)
% plot(t, p_c8)
% plot(t, p_c9)
% plot(t, p_c10)
plot(t, average_pressure)
plot(t_training_set, pressure_training_set)
set(gca,'FontSize',30)
xlabel('Time [s]')
ylabel('Pressure [Bar]')
legend('average','training set')
title('Pressure over Time E0 NL', 'FontSize', 50)

a   = 0.027;      %radius (equal to stroke/2)
B   = 0.068;      %diameter of piston
L   = 0.085;      %length of rod
V_c = 2.6148e-05; %deathvolume or clearance volume

s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

figure()
hold on
plot(Volume, average_pressure)
plot(Volume, pressure_training_set)
set(gca,'FontSize',30)
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
legend('average','training set')
title('E0 NL PV-diagram', 'FontSize', 50)
%% -------------------------HALF LOAD E0---------------------------------
clear all;close all;clc;

for i =0:9;
  
    opts = delimitedTextImportOptions("NumVariables", 3);
    
    rotation_speed = 47.08198676; %1 rotation = 0.0212395455s
    samples        = fix((1/rotation_speed)*100000*2);
    a              = samples * i;
    opts.DataLines = [1+a, samples+a];

    opts.Delimiter        = "\t";
    opts.VariableNames    = ["E0", "E1", "E2"];
    opts.VariableTypes    = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule    = "read";

    Test = readtable("D:\Documents\MATLAB\(4GB10) DBL Combustion engine\Experiment 2\raw_data_exp2\TestE0_half_load_3.txt", opts); 

    t          = Test.E0;
    puls_sens  = Test.E1;
    pres_sens  = Test.E2;
    
    x = -16; %correcting pressure sensor
    
    while i == 0
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c1 = pressure_relative + x;
        i = i+1
    end
    
    while i == 1
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c2 = pressure_relative + x;
        i = i+1
    end
    
    while i == 2
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c3 = pressure_relative + x;
        i = i+1
    end
    
    while i == 3
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c4 = pressure_relative + x;
        i = i+1
    end
    
    while i == 4
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c5 = pressure_relative + x;
        i = i+1
    end
    
    while i == 5
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c6 = pressure_relative + x;
        i = i+1
    end
    
    while i == 6
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c7 = pressure_relative + x;
        i = i+1
    end
    
    while i == 7
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c8 = pressure_relative + x;
        i = i+1
    end
    
    while i == 8
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c9 = pressure_relative + x;
        i = i+1
    end
    
    while i == 9
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c10 = pressure_relative + x;
        i = i+1
    end
    
    clear opts
end

average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

start_angle  = 210.852;
double_tooth = 83.13;
phi_speed    = rotation_speed*360*t; 
phi          = phi_speed-(start_angle+double_tooth);

%Training set half load
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines        = [6300, 10546];
opts.Delimiter        = "\t";
opts.VariableNames    = ["E0", "E1", "E2"];
opts.VariableTypes    = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule    = "read";

Training_set = readtable("D:\Documents\MATLAB\(4GB10) DBL Combustion engine\Test data\Training set\Half_load.txt", opts);

t_relative    = Training_set.E0;
puls_sens_ts  = Training_set.E2;
pres_sens_ts  = Training_set.E1;

pressure_relative_training_set = (pres_sens_ts -(0.115*5))/(0.00385*5*4); %formula from Canvas
pressure_training_set = pressure_relative_training_set*3.5-0.161;
t_training_set = t_relative + 0.320;

figure()
hold on 
% plot(t, p_c1)
% plot(t, p_c2)
% plot(t, p_c3)
% plot(t, p_c4)
% plot(t, p_c5)
% plot(t, p_c6)
% plot(t, p_c7)
% plot(t, p_c8)
% plot(t, p_c9)
% plot(t, p_c10)
plot(t, average_pressure)
plot(t_training_set, pressure_training_set)
set(gca,'FontSize',30)
xlabel('Time [s]')
ylabel('Pressure [Bar]')
legend('average', 'training set')
title('Pressure over Time E0 HL', 'FontSize', 50)

a   = 0.027;      %radius (equal to stroke/2)
B   = 0.068;      %diameter of piston
L   = 0.085;      %length of rod
V_c = 2.6148e-05; %deathvolume or clearance volume

s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

figure()
hold on
plot(Volume, average_pressure)
plot(Volume, pressure_training_set)
set(gca,'FontSize',30)
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
legend('average', 'training set')
title('E0 HL PV-diagram', 'FontSize', 50)
%% -------------------------FULL LOAD E0---------------------------------
clear all;close all;clc;

for i =0:9;
  
    opts = delimitedTextImportOptions("NumVariables", 3);
    
    rotation_speed = 45.17683503; %(1 rotation = 0.0221352381s
    samples        = fix((1/rotation_speed)*100000*2);
    a              = samples * i;
    opts.DataLines = [1+a, samples+a];

    opts.Delimiter        = "\t";
    opts.VariableNames    = ["E0", "E1", "E2"];
    opts.VariableTypes    = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule    = "read";

    Test = readtable("D:\Documents\MATLAB\(4GB10) DBL Combustion engine\Experiment 2\raw_data_exp2\TestE0_full_load_3.txt", opts); 

    t          = Test.E0;
    puls_sens  = Test.E1;
    pres_sens  = Test.E2;
    
    x = -24.9; %correcting pressure sensor
    
    while i == 0
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c1 = pressure_relative + x;
        i = i+1
    end
    
    while i == 1
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c2 = pressure_relative + x;
        i = i+1
    end
    
    while i == 2
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c3 = pressure_relative + x;
        i = i+1
    end
    
    while i == 3
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c4 = pressure_relative + x;
        i = i+1
    end
    
    while i == 4
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c5 = pressure_relative + x;
        i = i+1
    end
    
    while i == 5
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c6 = pressure_relative + x;
        i = i+1
    end
    
    while i == 6
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
        p_c7 = pressure_relative + x;
        i = i+1
    end
    
    while i == 7
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c8 = pressure_relative + x;
        i = i+1
    end
    
    while i == 8
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c9 = pressure_relative + x;
        i = i+1
    end
    
    while i == 9
        pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
        p_c10 = pressure_relative + x;
        i = i+1
    end
    
    clear opts
end

average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

start_angle  = 309.9854;
double_tooth = 83.13;
phi_speed    = rotation_speed*360*t; 
phi          = phi_speed-(start_angle+double_tooth);

%Training set Full load
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines        = [1700, 6126];
opts.Delimiter        = "\t";
opts.VariableNames    = ["E0", "E1", "E2"];
opts.VariableTypes    = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule    = "read";

Training_set = readtable("D:\Documents\MATLAB\(4GB10) DBL Combustion engine\Test data\Training set\Full_load.txt", opts);

t_relative    = Training_set.E0;
puls_sens_ts  = Training_set.E2;
pres_sens_ts  = Training_set.E1;

pressure_relative_training_set = (pres_sens_ts -(0.115*5))/(0.00385*5*4); %formula from Canvas
pressure_training_set = pressure_relative_training_set*2.9+0.1;
t_training_set = t_relative + 0.3814;

figure()
hold on 
% plot(t, p_c1)
% plot(t, p_c2)
% plot(t, p_c3)
% plot(t, p_c4)
% plot(t, p_c5)
% plot(t, p_c6)
% plot(t, p_c7)
% plot(t, p_c8)
% plot(t, p_c9)
% plot(t, p_c10)
plot(t, average_pressure)
plot(t_training_set, pressure_training_set)
set(gca,'FontSize',30)
xlabel('Time [s]')
ylabel('Pressure [Bar]')
legend('average', 'training set')
title('Pressure over Time E0 FL', 'FontSize', 50)

a   = 0.027;      %radius (equal to stroke/2)
B   = 0.068;      %diameter of piston
L   = 0.085;      %length of rod
V_c = 2.6148e-05; %deathvolume or clearance volume

s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

figure()
hold on
plot(Volume, average_pressure)
plot(Volume, pressure_training_set)
set(gca,'FontSize',30)
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
legend('average', 'training set')
title('E0 FL PV-diagram', 'FontSize', 50)