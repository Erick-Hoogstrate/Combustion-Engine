function [] = test(percentage, load)

if percentage == 0
    if load == 'no'
        %% -------------------------NO LOAD E0---------------------------------
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

            Test = readtable("raw_data_exp2\TestE0_no_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -1.4; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 72.1419;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)
        
    elseif load == 'half'
        %% -------------------------HALF LOAD E0---------------------------------
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

            Test = readtable("raw_data_exp2\TestE0_half_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -16; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 210.852;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)
        
    elseif load == 'full'
        %% -------------------------FULL LOAD E0---------------------------------
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

            Test = readtable("raw_data_exp2\TestE0_full_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -24.9; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 309.9854;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);
        
        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)
        
    end 
elseif percentage == 5
    if load == 'no'
        %% -------------------------NO LOAD E5---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 49.64010921; %1 rotation = 0.020145s  
            samples        = fix((1/rotation_speed)*100000*2); 
            a              = samples * i;
            opts.DataLines = [1+a, samples+a]; 

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE5_no_load_2.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = 0.1; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 63.6188;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    elseif load == 'half'
        %% -------------------------HALF LOAD E5---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 48.06486665; %1 rotation = 0.0208052174s
            samples        = fix((1/rotation_speed)*100000*2);
            a              = samples * i;
            opts.DataLines = [1+a, samples+a];

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE5_half_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -9.5; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 162.3054;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    elseif load == 'full'
        %% -------------------------FULL LOAD E5---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 45.40049719; %1 rotation = 0.0220261905s 
            samples        = fix((1/rotation_speed)*100000*2);
            a              = samples * i;
            opts.DataLines = [1+a, samples+a];

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE5_full_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -10.4; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 293.8683;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    end
elseif percentage == 10
    if load == 'no'
        %% -------------------------NO LOAD E10---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 49.73295566; %1 rotation = 0.0201073913s 
            samples        = fix((1/rotation_speed)*100000*2); 
            a              = samples * i;
            opts.DataLines = [1+a, samples+a]; 

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE10_no_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -0.4; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 356.1079;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    elseif load == 'half'
        %% -------------------------HALF LOAD E10---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 47.85883726; %1 rotation = 0.0208947826s 
            samples        = fix((1/rotation_speed)*100000*2);
            a              = samples * i;
            opts.DataLines = [1+a, samples+a];

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE10_half_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -6.7; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 265.674;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    elseif load == 'full'
        %% -------------------------FULL LOAD E10---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 45.88258429; %1 rotation = 0.0217947619s
            samples        = fix((1/rotation_speed)*100000*2);
            a              = samples * i;
            opts.DataLines = [1+a, samples+a];

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE10_full_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -14.2; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 345.2206;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    end
elseif percentage == 15
    if load == 'no'
        %% -------------------------NO LOAD E15---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 49.74814999; %1 rotation = 0.02010125s
            samples        = fix((1/rotation_speed)*100000*2);
            a              = samples * i;
            opts.DataLines = [1+a, samples+a];

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE15_no_load_1.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = 1.8+0.35; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 225.2637;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    elseif load == 'half'
        %% -------------------------HALF LOAD E15---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 48.15140473; %1 rotation = 0.0207678261s
            samples        = fix((1/rotation_speed)*100000*2);
            a              = samples * i;
            opts.DataLines = [1+a, samples+a];

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE15_half_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -8.9; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 204.0271;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)

    elseif load == 'full'
        %% -------------------------FULL LOAD E15---------------------------------
        for i =0:9;

            opts = delimitedTextImportOptions("NumVariables", 3);

            rotation_speed = 45.71713554; %1 rotation = 0.0218736364s 
            samples        = fix((1/rotation_speed)*100000*2);
            a              = samples * i;
            opts.DataLines = [1+a, samples+a];

            opts.Delimiter        = "\t";
            opts.VariableNames    = ["E0", "E1", "E2"];
            opts.VariableTypes    = ["double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule    = "read";

            Test = readtable("raw_data_exp2\TestE15_full_load_3.txt", opts); 

            t          = Test.E0;
            puls_sens  = Test.E1;
            pres_sens  = Test.E2;

            x = -16.9; %correcting pressure sensor

            while i == 0
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c1 = pressure_relative + x;
                i = i+1;
            end

            while i == 1
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c2 = pressure_relative + x;
                i = i+1;
            end

            while i == 2
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c3 = pressure_relative + x;
                i = i+1;
            end

            while i == 3
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c4 = pressure_relative + x;
                i = i+1;
            end

            while i == 4
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c5 = pressure_relative + x;
                i = i+1;
            end

            while i == 5
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c6 = pressure_relative + x;
                i = i+1;
            end

            while i == 6
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4);  
                p_c7 = pressure_relative + x;
                i = i+1;
            end

            while i == 7
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c8 = pressure_relative + x;
                i = i+1;
            end

            while i == 8
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c9 = pressure_relative + x;
                i = i+1;
            end

            while i == 9
                pressure_relative = (pres_sens -(0.115*5))/(0.00385*5*4); 
                p_c10 = pressure_relative + x;
                i = i+1;
            end

            clear opts
        end

        average_pressure = (p_c1 + p_c2 + p_c3 + p_c4 + p_c5 + p_c6 + p_c7 + p_c8 + p_c9 + p_c10)/10;

        start_angle  = 111.9155;
        double_tooth = 83.13;
        phi_speed    = rotation_speed*360*t; 
        phi          = phi_speed-(start_angle+double_tooth);

        a   = 0.027;      %radius (equal to stroke/2)
        B   = 0.068;      %diameter of piston
        L   = 0.085;      %length of rod
        V_c = 2.6148e-05; %deathvolume or clearance volume

        s      = a*cosd(phi) + (L.^2 + (a.^2)*(sind(phi)).^2).^(0.5); 
        Volume = (((pi)*B.^2)/(4))*(L + a - s) + V_c; 

        hold on
        plot(Volume, average_pressure,'LineWidth',6)
    end
    
else
    disp("Not included")
 
    
end %end if-elseif statement


end %end function