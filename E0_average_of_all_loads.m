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

    Test = readtable("C:\Users\20192420\Documents\Q3 Y2\data exp 2\TestE0_no_load_3.txt", opts); 

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

figure()
hold on 
plot(t, p_c1)
plot(t, p_c2)
plot(t, p_c3)
plot(t, p_c4) 
plot(t, p_c5) 
plot(t, p_c6)
plot(t, p_c7)
plot(t, p_c8)
plot(t, p_c9)
plot(t, p_c10)
plot(t, average_pressure)
xlabel('Time [s]')
ylabel('Pressure [Bar]')
legend('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10', 'average')
title('Pressure over Time E0 NL (10 cycles)')

a   = 0.027;      %radius (equal to stroke/2)
B   = 0.068;      %diameter of piston
L   = 0.085;      %length of rod
V_c = 2.6148e-05; %deathvolume or clearance volume

s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

figure()
plot(Volume, average_pressure)
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
title('E0 NL Average PV-diagram')
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

    Test = readtable("C:\Users\20192420\Documents\Q3 Y2\data exp 2\TestE0_half_load_3.txt", opts); 

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

figure()
hold on 
plot(t, p_c1)
plot(t, p_c2)
plot(t, p_c3)
plot(t, p_c4)
plot(t, p_c5)
plot(t, p_c6)
plot(t, p_c7)
plot(t, p_c8)
plot(t, p_c9)
plot(t, p_c10)
plot(t, average_pressure)
xlabel('Time [s]')
ylabel('Pressure [Bar]')
legend('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10', 'average')
title('Pressure over Time E0 HL (10 cycles)')

a   = 0.027;      %radius (equal to stroke/2)
B   = 0.068;      %diameter of piston
L   = 0.085;      %length of rod
V_c = 2.6148e-05; %deathvolume or clearance volume

s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

figure()
plot(Volume, average_pressure)
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
title('E0 HL Average PV-diagram')
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

    Test = readtable("C:\Users\20192420\Documents\Q3 Y2\data exp 2\TestE0_full_load_3.txt", opts); 

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

figure()
hold on 
plot(t, p_c1)
plot(t, p_c2)
plot(t, p_c3)
plot(t, p_c4)
plot(t, p_c5)
plot(t, p_c6)
plot(t, p_c7)
plot(t, p_c8)
plot(t, p_c9)
plot(t, p_c10)
plot(t, average_pressure)
xlabel('Time [s]')
ylabel('Pressure [Bar]')
legend('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10', 'average')
title('Pressure over Time E0 FL (10 cycles)')

a   = 0.027;      %radius (equal to stroke/2)
B   = 0.068;      %diameter of piston
L   = 0.085;      %length of rod
V_c = 2.6148e-05; %deathvolume or clearance volume

s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

figure()
plot(Volume, average_pressure)
xlabel('Volume [m^3]')
ylabel('Pressure [Bar]')
title('E0 FL Average PV-diagram')